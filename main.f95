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
    
    !Grandeurs caractériqtiques pour le calcul des nombres adimensionnels
    REAL(KIND = RKind) :: l_c, u_c
    

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
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: a, a_loc, a_opti
    !vecteur b de l'équation de poisson
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: b
    !Stockage de p sous forme vectorielle et residu de la méthode de jacobi
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: p_vec, p_vec_temp, residual, conjugate
    
    !structure qui contient les coordonnées en x et y d'un point du maillage
    TYPE COORDS
        REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: x, y
        INTEGER, DIMENSION(:, :), ALLOCATABLE :: borders
        REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: grad_x, grad_y
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
    
    !Permet de savoir de quelle manière appliquer des conditions initiales
    CHARACTER(LEN = 2) :: initial_conditions
    
    TYPE SETUP_TYPE
        REAL(KIND = RKind), DIMENSION(5) :: u, v
        REAL(KIND = RKind), DIMENSION(2, 5) :: poly
        REAL(KIND = RKind), DIMENSION(2) :: poly_ref, squares_ref
        REAL(KIND = RKind), DIMENSION(2, 3) :: squares
    END TYPE
    
    TYPE(SETUP_TYPE) :: setup
    
CONTAINS
    

    !maj à l'Etape 7, 2D, retiré les variables plus utilisées et ajout du choix de scheme pour les termes convectifs
    SUBROUTINE read_input_file(name)
    IMPLICIT NONE
        
        CHARACTER(LEN = StrLen), INTENT(IN) :: name
        INTEGER :: i
        
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
        READ(10, *)
        READ(10, *) l_c
        READ(10, *)
        READ(10, *) u_c
        READ(10, *)
        viscosity = u_c*l_c/re
        
        !Scheme choisi pour termes convectifs
        READ(10, *)
        READ(10, *) scheme

        READ(10, *)
        READ(10, *)
        READ(10, *)
        READ(10, *)
        READ(10, *) initial_conditions
        
        IF (initial_conditions == 'CC') THEN
            
            READ(10, *)
            READ(10, *)
            READ(10, *) setup%u(4)
            READ(10, *) setup%u(1), setup%u(2)
            READ(10, *) setup%u(3)
            READ(10, *)
            READ(10, *) setup%v(4)
            READ(10, *) setup%v(1), setup%v(2)
            READ(10, *) setup%v(3)
            READ(10, *)
            READ(10, *) setup%u(5), setup%v(5)
            
            READ(10, *)
            
            READ(10, *)
            READ(10, *)
            READ(10, *) setup%poly(1, :)
            READ(10, *) setup%poly(2, :)
            READ(10, *)
            READ(10, *) setup%poly_ref(:)
            
            READ(10, *)
            
            READ(10, *)
            READ(10, *) setup%squares(1, :)
            READ(10, *) setup%squares(2, :)
            READ(10, *)
            READ(10, *) setup%squares_ref(:)
            
        ELSE
            
            setup%u(1:3) = 0_RKind
            setup%u(4) = 1_RKind
            setup%u(5) = 0.0_RKind
            setup%v(:) = 0_RKind
            setup%poly(:, :4) = 0_RKind
            setup%poly(:, 5) = -1_RKind
            setup%poly_ref(:) = 0.0_RKind
            setup%poly_ref(:) = 0.0_RKind
            setup%squares(2, :2) = 0_RKind
            setup%squares(:, 3) = -1_RKind
            setup%squares_ref(:) = 0.0_RKind
            setup%squares_ref(:) = 0.0_RKind
            
        END IF
        
        
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
            '   seconds", I=', n_y, ', J=', n_x, ', DATAPACKING=POINT'
        
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
        name = 'debug/a_opti_' // TRIM(name) // '.dat'
        
        !Ouverture du fichier a ecrire
        OPEN(11, FILE = name)
        
        DO i = 1, n_x*n_y
            WRITE(11, '(A, I3, A, 5(F7.1))') 'i = ', i, ' | a_opti(i) = ', a_opti(i, :)
        END DO
        
        !Fermeture du fichier
        CLOSE(11)
        
    END SUBROUTINE Debug
    
        
    
    
    !maj à l'Etape 7, 2D
    !subroutine de creation du maillage spatial, contient les coordonnées exactes de chaque pt, en 2D
    SUBROUTINE create_space_grid()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j, k
        INTEGER, DIMENSION(:, :), ALLOCATABLE :: borders_grid
        REAL(KIND = RKIND) :: x, y, poly1, poly2, squares1, squares2
        
        ALLOCATE(space_grid%x(n_x))
        ALLOCATE(space_grid%y(n_y))
        ALLOCATE(space_grid%borders(n_x, n_y))
        ALLOCATE(borders_grid(-1:n_x+2, -1:n_y+2))
        ALLOCATE(space_grid%grad_x(0:n_x+1, 0:n_y+1))
        ALLOCATE(space_grid%grad_y(0:n_x+1, 0:n_y+1))
        
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
        
        borders_grid(:, :) = 0
        borders_grid(-1:0, :) = 1
        borders_grid(n_x+1:n_x+2, :) = 1
        borders_grid(:, -1:0) = 1
        borders_grid(:, n_y+1:n_y+2) = 1
        
        ! OPEN(11, FILE = 'debug/borders_grid.txt')
        
        ! DO i = -1, n_x+2
        !     WRITE(11, *) borders_grid(i, :)
        ! END DO
        
        ! CLOSE(11)
        
        space_grid%grad_x(:, :) = 0.0_RKind
        space_grid%grad_y(:, :) = 0.0_RKind
        
        space_grid%grad_y(:, 0) = 1.0_RKind
        space_grid%grad_y(:, n_y+1) = -1.0_RKind
        space_grid%grad_x(0, :) = 1.0_RKind
        space_grid%grad_x(n_x+1, :) = -1.0_RKind
        
        DO j = 1, n_y
            y = space_grid%y(j)
            DO i = 1, n_x
                x = space_grid%x(i)
                poly1 = (x - setup%poly_ref(1))*(setup%poly(1, 1)*(x - setup%poly_ref(1)) + setup%poly(1, 2)) &
                + (y - setup%poly_ref(2))*(setup%poly(1, 3)*(y - setup%poly_ref(2)) + setup%poly(1, 4))
                poly1 = setup%poly(1, 5) - poly1
                
                poly2 = (x - setup%poly_ref(1))*(setup%poly(2, 1)*(x - setup%poly_ref(1)) + setup%poly(2, 2)) &
                + (y - setup%poly_ref(2))*(setup%poly(2, 3)*(y - setup%poly_ref(2)) + setup%poly(2, 4))
                poly2 = setup%poly(2, 5) - poly2
                
                squares1 = ABS(x - setup%squares_ref(1))*setup%squares(1, 1) + ABS(y - setup%squares_ref(2))*setup%squares(1, 2)
                squares1 = setup%squares(1, 3) - squares1
                
                squares2 = ABS(x - setup%squares_ref(1))*setup%squares(2, 1) + ABS(y - setup%squares_ref(2))*setup%squares(2, 2)
                squares2 = setup%squares(2, 3) - squares2
                
                IF ((poly1 >= 0) .AND. (poly2 >= 0) .AND. (squares1 >= 0) .AND. (squares2 >= 0)) THEN
                
                    borders_grid(i, j) = 1
                    
                    IF ((poly1 <= poly2) .AND. (poly1 <= squares1) .AND. (poly1 <= squares2)) THEN
                        
                        space_grid%grad_x(i, j) = 2*setup%poly(1, 1)*(x - setup%poly_ref(1)) + setup%poly(1, 2)
                        space_grid%grad_y(i, j) = 2*setup%poly(1, 3)*(y - setup%poly_ref(2)) + setup%poly(1, 4)
                        
                    ELSE IF ((poly2 <= squares1) .AND. (poly2 <= squares2)) THEN
                        
                        space_grid%grad_x(i, j) = 2*setup%poly(2, 1)*(x - setup%poly_ref(1)) + setup%poly(2, 2)
                        space_grid%grad_y(i, j) = 2*setup%poly(2, 3)*(y - setup%poly_ref(2)) + setup%poly(2, 4)
                        
                    ELSE IF (squares1 <= squares2) THEN
                        
                        space_grid%grad_x(i, j) = setup%squares(1, 1)*(x - setup%squares_ref(1))/ABS((x - setup%squares_ref(1)))
                        space_grid%grad_y(i, j) = setup%squares(1, 2)*(y - setup%squares_ref(2))/ABS((y - setup%squares_ref(2)))
                        
                    ELSE
                        
                        space_grid%grad_x(i, j) = setup%squares(2, 1)*(x - setup%squares_ref(1))/ABS((x - setup%squares_ref(1)))
                        space_grid%grad_y(i, j) = setup%squares(2, 2)*(y - setup%squares_ref(2))/ABS((y - setup%squares_ref(2)))
                        
                    END IF
                    
                END IF
                
            END DO
        END DO
        
        ! OPEN(11, FILE = 'debug/borders_grid.txt')
        
        ! DO i = 1, n_x
        !     WRITE(11, *) borders_grid(i, :)
        ! END DO
        
        ! CLOSE(11)
        
        DO j = 1, n_y
            space_grid%borders(:, j) = borders_grid(0:n_x-1, j)*1 + borders_grid(-1:n_x-1, j)*16 &
            + borders_grid(2:n_x+1, j)*2 + borders_grid(3:n_x+2, j)*32 &
            + borders_grid(1:n_x, j-1)*4 + borders_grid(1:n_x, j-2)*64 &
            + borders_grid(1:n_x, j+1)*8 + borders_grid(1:n_x, j+2)*128
        END DO
        
        DO j = 1, n_y
            DO i = 1, n_x
                IF (borders_grid(i, j) == 1) THEN
                    space_grid%borders(i, j) = -1
                END IF
            END DO
        END DO
        
        OPEN(11, FILE = 'debug/borders.txt')
        
        DO i = 1, n_x
            WRITE(11, *) space_grid%borders(i, :)
        END DO
        
        CLOSE(11)
                
        !j\i     -2  -1  0   +1  +2
            !+2  . | . |128| . | . 
            !+1  . | . | 8 | . | . 
            ! 0  16| 1 | x | 2 | 32
            !-1  . | . | 4 | . | . 
            !-2  . | . | 64| . | . 
        
        DEALLOCATE(borders_grid)
        
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
        u(:, :) = setup%u(5)
        v(:, :) = setup%v(5)
        p(:, :) = 0_RKind
        
        u(1, :) = setup%u(1)
        u(n_x, :) = setup%u(2)
        u(:, 1) = setup%u(3)
        u(:, n_y) = setup%u(4)
        
        v(1, :) = setup%v(1)
        v(n_x, :) = setup%v(2)
        v(:, 1) = setup%v(3)
        v(:, n_y) = setup%v(4)
        
        DO j = 1, n_y
            DO i = 1, n_x
                IF (space_grid%borders(i, j) < 0) THEN
                    u(i, j) = 0_RKind
                    v(i, j) = 0_RKind
                    p(i, j) = 0_RKind
                END IF
            END DO
        END DO
        
        
        
    END SUBROUTINE init_solution
    
    
    
    !maj à l'étape 4, 2D
    !Création et remplissage de la matrice A et allocation des tableaux nécessaires à Jacobi
    SUBROUTINE init_a()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j, k, k_max
        REAL(KIND = RKIND) :: inv_x_2, inv_y_2

        !calcul de la dimension du vecteur solution
        k_max = n_x*n_y
        
        ALLOCATE(a(k_max, k_max))
        ALLOCATE(a_loc(k_max, k_max))
        ALLOCATE(a_opti(k_max, 5))
        ALLOCATE(p_vec(k_max))
        ALLOCATE(p_vec_temp(k_max))
        ALLOCATE(residual(k_max))
        ALLOCATE(b(k_max))
        ALLOCATE(conjugate(k_max))
        
        
        
        !Remplissage de la matrice A
        !initialisation
        a(:, :) = 0
        
        !calcul des coefficients
        inv_x_2 = 1_RKind/(dx**2.0_RKind)
        inv_y_2 = 1_RKind/(dy**2.0_RKind)
        
        !j\i     -2  -1  0   +1  +2
            !+2  . | . |128| . | . 
            !+1  . | . | 8 | . | . 
            ! 0  16| 1 | x | 2 | 32
            !-1  . | . | 4 | . | . 
            !-2  . | . | 64| . | .
        !
        
        DO j = 1, n_y
            DO i = 1, n_x

                k = (j - 1)*n_x + i
                
                IF (space_grid%borders(i,j) < 0) THEN
                
                    a(k, k) = 0_RKind
                    
                !Conditions limites de dérivée nulle
                ELSE IF (MOD(space_grid%borders(i,j), 2) == 1) THEN
                    
                    a(k, k) = -ABS(space_grid%grad_x(i-1, j))/dx - ABS(space_grid%grad_y(i-1, j))/dy
                    a(k, k + 1) = ABS(space_grid%grad_x(i-1, j))/dx
                    IF (space_grid%grad_y(i-1, j) > 0) THEN
                        a(k, k + n_x) = ABS(space_grid%grad_y(i-1, j))/dy
                    ELSE
                        a(k, k - n_x) = ABS(space_grid%grad_y(i-1, j))/dy
                    END IF

                ELSE IF (MOD(space_grid%borders(i,j), 4)/2 == 1) THEN
                    
                    a(k, k) = -ABS(space_grid%grad_x(i+1, j))/dx - ABS(space_grid%grad_y(i+1, j))/dy
                    a(k, k - 1) = ABS(space_grid%grad_x(i+1, j))/dx
                    IF (space_grid%grad_y(i+1, j) > 0) THEN
                        a(k, k + n_x) = ABS(space_grid%grad_y(i+1, j))/dy
                    ELSE
                        a(k, k - n_x) = ABS(space_grid%grad_y(i+1, j))/dy
                    END IF
                    
                ELSE IF (MOD(space_grid%borders(i,j), 8)/4 == 1) THEN
                
                    a(k, k) = -ABS(space_grid%grad_x(i, j-1))/dx - ABS(space_grid%grad_y(i, j-1))/dy
                    a(k, k + n_x) = ABS(space_grid%grad_y(i, j-1))/dy
                    IF (space_grid%grad_x(i, j-1) > 0) THEN
                        a(k, k + 1) = ABS(space_grid%grad_x(i, j-1))/dx
                    ELSE
                        a(k, k - 1) = ABS(space_grid%grad_x(i, j-1))/dx
                    END IF
                    
                ELSE IF (MOD(space_grid%borders(i,j), 16)/8 == 1) THEN
                    
                    a(k, k) = -ABS(space_grid%grad_x(i, j+1))/dx - ABS(space_grid%grad_y(i, j+1))/dy
                    a(k, k - n_x) = ABS(space_grid%grad_y(i, j+1))/dy
                    IF (space_grid%grad_x(i, j+1) > 0) THEN
                        a(k, k + 1) = ABS(space_grid%grad_x(i, j+1))/dx
                    ELSE
                        a(k, k - 1) = ABS(space_grid%grad_x(i, j+1))/dx
                    END IF
                    
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
        
        a_opti(:, :) = 0.0_RKind
        DO j = 1, n_y
            DO i = 1, n_x

                k = (j - 1)*n_x + i
                
                a_opti(k, 3) = a(k, k)
                IF (i > 1) THEN
                    a_opti(k, 2) = a(k, k-1)
                END IF
                IF (i < n_x) THEN
                    a_opti(k, 4) = a(k, k+1)
                END IF
                IF (j > 0) THEN
                    a_opti(k, 1) = a(k, k-n_x)
                END IF
                IF (j < n_y) THEN
                    a_opti(k, 5) = a(k, k+n_x)
                END IF
                
                ! IF (SUM(a_opti(k, :)) /= 0) THEN
                !     PRINT*, 'i = ', i, 'j = ', j, SUM(a_opti(k, :))
                ! END IF
                
            END DO
        END DO
        
    END SUBROUTINE init_a
    
    
    
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
        
        
        
        !création du maillage spatial
        CALL create_space_grid()
    
        !initalisation de la solution grâce aux conditions initiales
        CALL init_solution()
        
        CALL init_a()
        
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
                    !PRINT*, dt_temp, 'u'
                END IF

                dt_temp = cfl*dy/ABS(v(i,j))     !calcul du minimum potentiel sur v

                IF (dt_temp < dt_min) THEN
                    dt_min = dt_temp
                    !PRINT*, dt_temp, 'v'
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
                IF (MOD(space_grid%borders(i, j), 16) /= 0) THEN

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
    
    
    !Integrale de la pression
    SUBROUTINE pressure_integral_correction(integral)
    
    IMPLICIT NONE
        
        REAL(KIND = RKind) :: integral
        INTEGER(KIND = IKIND) :: i, j
        
        integral = 0_RKIND
        DO i = 1, n_x
            DO j = 1, n_y
                IF (((i == 1) .OR. (i == n_x)) .NEQV. ((j == 1) .OR. (j == n_y))) THEN
                    integral = integral + p_vec((j-1)*n_x+i)*dx*dy*0.25_RKind
                ELSE IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN
                    integral = integral + p_vec((j-1)*n_x+i)*dx*dy*0.5_RKind
                ELSE
                    integral = integral + p_vec((j-1)*n_x+i)*dx*dy
                END IF
            END DO
        END DO
        integral = integral/(l_x*l_y)
        p_vec(:) = p_vec(:) - integral
        
    END SUBROUTINE pressure_integral_correction
    
    
    
    !Methode de Jacobi (resolution iterative de systeme lineaire)
    SUBROUTINE jacobi_method()
    
    IMPLICIT NONE
        
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001     !point d'arret de jacobi
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000       !Arret forcé de jacobi
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        
        !CALL CPU_TIME(time1)
        
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        p_vec_temp(:) = p_vec(:)
        
        k_max = n_x*n_y
        
        residual(:) = - b(:)
        
        residual(1:k_max) = residual(1:k_max) + a_opti(1:k_max, 3)*p_vec(1:k_max)
        residual(2:k_max) = residual(2:k_max) + a_opti(2:k_max, 2)*p_vec(1:k_max-1)
        residual(1:k_max-1) = residual(1:k_max-1) + a_opti(1:k_max-1, 4)*p_vec(2:k_max)
        residual(n_x+1:k_max) = residual(n_x+1:k_max) + a_opti(n_x+1:k_max, 1)*p_vec(1:k_max-n_x)
        residual(1:k_max-n_x) = residual(1:k_max-n_x) + a_opti(1:k_max-n_x, 5)*p_vec(1+n_x:k_max)
        
        iteration = 0
        
        CALL pressure_integral_correction(integral)
        
        lower_norm = 1
        upper_norm = 1
        DO WHILE (upper_norm/lower_norm > RTol)
            
            p_vec_temp(:) = p_vec(:)
            DO i = 1, k_max
                IF (ABS(a_opti(i, 3)) > 1E-15_RKind) THEN
                    p_vec(i) = p_vec_temp(i) - residual(i)/a_opti(i, 3)
                END IF
            END DO
            
            CALL pressure_integral_correction(integral)

            CALL norm_2(p_vec, lower_norm)
            residual(:) = (p_vec(:)-p_vec_temp(:))
            CALL norm_2(residual, upper_norm)
            
            residual(:) = - b(:)
            
            residual(1:k_max) = residual(1:k_max) + a_opti(1:k_max, 3)*p_vec(1:k_max)
            residual(2:k_max) = residual(2:k_max) + a_opti(2:k_max, 2)*p_vec(1:k_max-1)
            residual(1:k_max-1) = residual(1:k_max-1) + a_opti(1:k_max-1, 4)*p_vec(2:k_max)
            residual(n_x+1:k_max) = residual(n_x+1:k_max) + a_opti(n_x+1:k_max, 1)*p_vec(1:k_max-n_x)
            residual(1:k_max-n_x) = residual(1:k_max-n_x) + a_opti(1:k_max-n_x, 5)*p_vec(1+n_x:k_max)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max de Jacobi atteint (', IterationMax, ' )'
                STOP
            END IF
            
            ! IF (MOD(iteration, 100) == 0) THEN
            !     PRINT*, 'iteration = ', iteration, ' | convergence = ', upper_norm/lower_norm
            ! END IF
            
        END DO
        
        !CALL CPU_TIME(time2)
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        !PRINT*, 'Jacobi for a grid size of ', n_x, ' : ', time2 - time1, ' seconds (', iteration, ' iterations)'
        PRINT*, 'Jacobi :', iteration, ', iterations | integrale(p) = ', integral
        
    END SUBROUTINE jacobi_method
    
    
    
    !Methode de Gauss-Siedel (resolution iterative de systeme lineaire)
    SUBROUTINE gauss_siedel_method()
    
    IMPLICIT NONE
        
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001     !point d'arret de jacobi
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000       !Arret forcé de jacobi
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        
        !CALL CPU_TIME(time1)
        
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        p_vec_temp(:) = p_vec(:)
        
        k_max = n_x*n_y
        
        iteration = 0
        
        CALL pressure_integral_correction(integral)
        
        
        lower_norm = 1
        upper_norm = 1
        DO WHILE (upper_norm/lower_norm > RTol)
            
            p_vec_temp(:) = p_vec(:)
            DO i = 1, k_max
                IF (ABS(a_opti(i, 3)) > 1E-15_RKind) THEN
                    p_vec(i) = b(i)
                    IF (i > 1) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 2)*p_vec(i-1)
                    END IF
                    IF (i > n_x) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 1)*p_vec(i-n_x)
                    END IF
                    IF (i < k_max) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 4)*p_vec(i+1)
                    END IF
                    IF (i <= k_max-n_x) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 5)*p_vec(i+n_x)
                    END IF
                    p_vec(i) = p_vec(i)/a_opti(i, 3)
                END IF
            END DO
            
            CALL pressure_integral_correction(integral)
            
            CALL norm_2(p_vec, lower_norm)
            residual(:) = (p_vec(:)-p_vec_temp(:))
            CALL norm_2(residual, upper_norm)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max de Gauss-Siedel atteint (', IterationMax, ' )'
                STOP
            END IF
            
            ! IF (MOD(iteration, 100) == 0) THEN
            !     PRINT*, 'iteration = ', iteration, ' | convergence = ', upper_norm/lower_norm
            ! END IF
            
        END DO
        
        !CALL CPU_TIME(time2)
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        !PRINT*, 'Jacobi for a grid size of ', n_x, ' : ', time2 - time1, ' seconds (', iteration, ' iterations)'
        PRINT*, 'Gauss-Siedel :', iteration, ', iterations | integrale(p) = ', integral
        
    END SUBROUTINE gauss_siedel_method
    
    
    
    !Methode de Jacobi (resolution iterative de systeme lineaire)
    SUBROUTINE successive_over_relaxation_method()
    
    IMPLICIT NONE
        
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001     !point d'arret de jacobi
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000       !Arret forcé de jacobi
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        REAL(KIND = RKind) :: sor_coeff
        
        sor_coeff = 1.4
        
        !CALL CPU_TIME(time1)
        
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        p_vec_temp(:) = p_vec(:)
        
        k_max = n_x*n_y
        
        iteration = 0
        
        CALL pressure_integral_correction(integral)
        
        
        lower_norm = 1
        upper_norm = 1
        DO WHILE (upper_norm/lower_norm > RTol)
            
            p_vec_temp(:) = p_vec(:)
            DO i = 1, k_max
                IF (ABS(a_opti(i, 3)) > 1E-15_RKind) THEN
                    p_vec(i) = b(i)
                    IF (i > 1) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 2)*p_vec(i-1)
                    END IF
                    IF (i > n_x) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 1)*p_vec(i-n_x)
                    END IF
                    IF (i < k_max) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 4)*p_vec(i+1)
                    END IF
                    IF (i <= k_max-n_x) THEN
                        p_vec(i) = p_vec(i) - a_opti(i, 5)*p_vec(i+n_x)
                    END IF
                    p_vec(i) = (1_RKind - sor_coeff)*p_vec_temp(i) + sor_coeff*p_vec(i)/a_opti(i, 3)
                END IF
            END DO
            
            CALL pressure_integral_correction(integral)
            
            CALL norm_2(p_vec, lower_norm)
            residual(:) = (p_vec(:)-p_vec_temp(:))
            CALL norm_2(residual, upper_norm)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max de Surrelaxation atteint (', IterationMax, ' )'
                STOP
            END IF
            
            ! IF (MOD(iteration, 100) == 0) THEN
            !     PRINT*, 'iteration = ', iteration, ' | convergence = ', upper_norm/lower_norm
            ! END IF
            
        END DO
        
        !CALL CPU_TIME(time2)
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        !PRINT*, 'Jacobi for a grid size of ', n_x, ' : ', time2 - time1, ' seconds (', iteration, ' iterations)'
        PRINT*, 'Surrelaxation :', iteration, ', iterations | integrale(p) = ', integral
        
    END SUBROUTINE successive_over_relaxation_method
    
    
    
    !Methode de Jacobi (resolution iterative de systeme lineaire)
    SUBROUTINE conjugate_gradient_method()
    
    IMPLICIT NONE
        
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001     !point d'arret de jacobi
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000       !Arret forcé de jacobi
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral, r_norm0, r_norm
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        REAL(KIND = RKind) :: alpha, beta
        
        !CALL CPU_TIME(time1)
        
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        p_vec_temp(:) = p_vec(:)
        
        k_max = n_x*n_y
        
        
        iteration = 0
        
        CALL pressure_integral_correction(integral)
        
        residual(:) = b(:) - MATMUL(a, p_vec)
        conjugate(:) = residual(:)
        
        CALL norm_2(residual, r_norm0)
        r_norm = r_norm0
        
        lower_norm = 1
        upper_norm = 1
        !DO WHILE (upper_norm/lower_norm > RTol)
        DO WHILE (r_norm > r_norm0*RTol)
            
            p_vec_temp(:) = p_vec(:)
            
            ! alpha = DOT_PRODUCT(residual, residual)/DOT_PRODUCT(conjugate, MATMUL(a, conjugate))
            ! p_vec(:) = p_vec_temp(:) + alpha*conjugate(:)
            ! beta = DOT_PRODUCT(residual, residual)
            ! residual(:) = residual(:) - alpha*MATMUL(a, conjugate)
            ! !residual(:) = -MATMUL(a, p_vec) + b(:)
            
            ! beta = DOT_PRODUCT(residual, residual)/beta
            
            ! conjugate(:) = residual(:) + beta*conjugate(:)
            
            alpha = DOT_PRODUCT(residual, residual)/DOT_PRODUCT(conjugate, MATMUL(a, conjugate))
            p_vec(:) = p_vec_temp(:) + alpha*conjugate(:)
            beta = DOT_PRODUCT(residual, residual)
            residual(:) = residual(:) - alpha*MATMUL(a, conjugate)
            !residual(:) = -MATMUL(a, p_vec) + b(:)
            
            beta = DOT_PRODUCT(residual, residual)/beta
            
            conjugate(:) = residual(:) + beta*conjugate(:)
            
            
            
            
            CALL pressure_integral_correction(integral)
            
            ! CALL norm_2(p_vec, lower_norm)
            ! p_vec_temp(:) = p_vec(:)-p_vec_temp(:)
            ! CALL norm_2(p_vec_temp, upper_norm)
            
            CALL norm_2(residual, r_norm)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max du Gradient conjugue atteint (', IterationMax, ' )'
                STOP
            END IF
            
            IF (MOD(iteration, 100) == 0) THEN
                PRINT*, 'iteration = ', iteration, ' | convergence = ', r_norm/r_norm0
            END IF
            
        END DO
        
        !CALL CPU_TIME(time2)
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        !PRINT*, 'Jacobi for a grid size of ', n_x, ' : ', time2 - time1, ' seconds (', iteration, ' iterations)'
        PRINT*, 'Gradient conjugue :', iteration, ', iterations | integrale(p) = ', integral
        
        
        
    END SUBROUTINE conjugate_gradient_method
    
    
    
    SUBROUTINE cd_scheme_2(i, j)
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind), INTENT(IN) :: i, j
        
        u_temp(i,j) = u(i,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind)) &
        + u(i-1,j)*dt*(viscosity/dx**2_RKind + u(i,j)/(2_RKind*dx)) &
        + u(i,j-1)*dt*(viscosity/dy**2_RKind + v(i,j)/(2_RKind*dy)) &
        + u(i+1,j)*dt*(viscosity/dx**2_RKind - u(i,j)/(2_RKind*dx)) &
        + u(i,j+1)*dt*(viscosity/dy**2_RKind - v(i,j)/(2_RKind*dy))
        
        !Calcul de v
        v_temp(i,j) = v(i,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind)) &
        + v(i-1,j)*dt*(viscosity/dx**2_RKind + u(i,j)/(2_RKind*dx)) &
        + v(i,j-1)*dt*(viscosity/dy**2_RKind + v(i,j)/(2_RKind*dy)) &
        + v(i+1,j)*dt*(viscosity/dx**2_RKind - u(i,j)/(2_RKind*dx)) &
        + v(i,j+1)*dt*(viscosity/dy**2_RKind - v(i,j)/(2_RKind*dy))
        
    END SUBROUTINE cd_scheme_2
    
    
    
    !calcul du profil de vitesse pour une itération temporelle
    SUBROUTINE speed_guess()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j !compteur
        INTEGER(KIND = IKIND) :: upwind_x, upwind_y     !permet de déterminer la direction du upwind (1 si backward, 0 sinon)
        
        
        !Conditions limites
        u_temp(1, :) = setup%u(1)
        u_temp(n_x, :) = setup%u(2)
        u_temp(:, 1) = setup%u(3)
        u_temp(:, n_y) = setup%u(4)
        
        v_temp(1, :) = setup%v(1)
        v_temp(n_x, :) = setup%v(2)
        v_temp(:, 1) = setup%v(3)
        v_temp(:, n_y) = setup%v(4)
        
        DO j = 1, n_y
            DO i = 1, n_x
                IF (space_grid%borders(i, j) < 0) THEN
                    u(i, j) = 0_RKind
                    v(i, j) = 0_RKind
                    p(i, j) = 0_RKind
                END IF
            END DO
        END DO
        
        
        
        !Calcul avec scheme régressif upwind d'ordre 1
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
        
        !Calcul avec scheme centré d'ordre 4
        ELSE IF (scheme == 'CD4') THEN
            
            !calcul du n+1
            
            
            DO i = 2, n_x-1
                CALL cd_scheme_2(i, 2)
                CALL cd_scheme_2(i, n_y - 1)
            END DO
            DO j = 3, n_y-2
                CALL cd_scheme_2(2, j)
                CALL cd_scheme_2(n_x-1, j)
            END DO
            
            DO j = 3, n_y-2
                    
                !Calcul de u
                u_temp(3:n_x-2,j) = u(3:n_x-2,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind)) &
                + u(2:n_x-3,j)*dt*(viscosity/dx**2_RKind + 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + u(3:n_x-2,j-1)*dt*(viscosity/dy**2_RKind + 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                + u(4:n_x-1,j)*dt*(viscosity/dx**2_RKind - 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + u(3:n_x-2,j+1)*dt*(viscosity/dy**2_RKind - 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                - u(1:n_x-4,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                - u(3:n_x-2,j-2)*dt*v(3:n_x-2,j)/(12_RKind*dx) &
                + u(5:n_x-1,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                + u(3:n_x-2,j+2)*dt*v(3:n_x-2,j)/(12_RKind*dx)
                
                !Calcul de v
                v_temp(3:n_x-2,j) = v(3:n_x-2,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind)) &
                + v(2:n_x-3,j)*dt*(viscosity/dx**2_RKind + 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + v(3:n_x-2,j-1)*dt*(viscosity/dy**2_RKind + 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                + v(4:n_x-1,j)*dt*(viscosity/dx**2_RKind - 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + v(3:n_x-2,j+1)*dt*(viscosity/dy**2_RKind - 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                - v(1:n_x-4,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                - v(3:n_x-2,j-2)*dt*v(3:n_x-2,j)/(12_RKind*dx) &
                + v(5:n_x-1,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                + v(3:n_x-2,j+2)*dt*v(3:n_x-2,j)/(12_RKind*dx)
                
            END DO
            
            !Schéma centré d'ordre 2 pour les bords
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    IF (space_grid%borders(i, j)/16 >= 0) THEN
                        CALL cd_scheme_2(i, j)
                    END IF
                END DO
            END DO
            
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    IF (MOD(space_grid%borders(i, j), 16) /= 0) THEN
                        u_temp(i, j) = 0_RKIND
                        v_temp(i, j) = 0_RKIND
                    END IF
                END DO
            END DO
            
        !Calcul avec scheme centré d'ordre 4
        ELSE IF (scheme == 'CD2') THEN
            
            !Schéma centré d'ordre 2 pour les bords
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    CALL cd_scheme_2(i, j)
                END DO
            END DO
            
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    IF (MOD(space_grid%borders(i, j), 16) /= 0) THEN
                        u_temp(i, j) = 0_RKIND
                        v_temp(i, j) = 0_RKIND
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
        !CALL gauss_siedel_method()
        !CALL successive_over_relaxation_method()
        !CALL conjugate_gradient_method()

    END SUBROUTINE compute_pressure
    
    
    SUBROUTINE adjust_speed()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j
        
        !Conditions limites
        u(1, :) = setup%u(1)
        u(n_x, :) = setup%u(2)
        u(:, 1) = setup%u(3)
        u(:, n_y) = setup%u(4)
        
        v(1, :) = setup%v(1)
        v(n_x, :) = setup%v(2)
        v(:, 1) = setup%v(3)
        v(:, n_y) = setup%v(4)
        
        DO j = 1, n_y
            DO i = 1, n_x
                IF (space_grid%borders(i, j) < 0) THEN
                    u(i, j) = 0_RKind
                    v(i, j) = 0_RKind
                    p(i, j) = 0_RKind
                END IF
            END DO
        END DO
        
        DO j = 2, n_y-1
            u(2:n_x-1, j) = u_temp(2:n_x-1, j) - ((p(3:n_x, j) - p(1:n_x-2, j))/(2.0_RKind*dx))*dt/density
            v(2:n_x-1, j) = v_temp(2:n_x-1, j) - ((p(2:n_x-1, j+1) - p(2:n_x-1, j-1))/(2.0_RKind*dy))*dt/density
        END DO
        
        DO j = 2, n_y-1
            DO i = 2, n_x-1
                IF (MOD(space_grid%borders(i, j), 16) /= 0) THEN
                    u(i, j) = 0_RKIND
                    v(i, j) = 0_RKIND
                END IF
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
    
    REAL(KIND = RKIND) :: time1, time2
    
    CALL CPU_TIME(time1)
    
    !récupération des données du problème
    name = 'input.dat'
    CALL read_input_file(name)
    
    CALL initialisation()
    
    CALL debug(0)
    
    CALL resolution_loop()
    
    CALL CPU_TIME(time2)
    
    PRINT*, 'temps pour la résolution : ', time2-time1
    
    DEALLOCATE(space_grid%x)
    DEALLOCATE(space_grid%y)
    DEALLOCATE(space_grid%borders)
    DEALLOCATE(u)
    DEALLOCATE(v)
    DEALLOCATE(a)
    DEALLOCATE(a_loc)
    DEALLOCATE(a_opti)
    DEALLOCATE(b)
    DEALLOCATE(p)
    DEALLOCATE(p_vec)
    DEALLOCATE(p_vec_temp)
    DEALLOCATE(residual)
    DEALLOCATE(space_grid%grad_x)
    DEALLOCATE(space_grid%grad_y)
    DEALLOCATE(conjugate)
    
END PROGRAM main