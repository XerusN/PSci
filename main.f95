MODULE global
!$ use OMP_LIB
IMPLICIT NONE

    !definition des tailles des reels(RKind) et des entiers(IKind)
    INTEGER, PARAMETER :: RKind = SELECTED_REAL_KIND(15,200)
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
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: a_opti
    !vecteur b de l'équation de poisson
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: b
    !Stockage de p sous forme vectorielle et residu de la méthode de jacobi
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: p_vec, p_vec_temp, residual, conjugate, a_mul_conj, a_mul_residual
    
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
    
    !Permet de définir les conditions initiales
    TYPE SETUP_TYPE
        REAL(KIND = RKind), DIMENSION(5) :: u, v
        REAL(KIND = RKind), DIMENSION(2, 5) :: poly
        REAL(KIND = RKind), DIMENSION(2) :: poly_ref, squares_ref
        REAL(KIND = RKind), DIMENSION(2, 3) :: squares
    END TYPE
    
    TYPE(SETUP_TYPE) :: setup
    
    !Variables pour le benchmark
    REAL(KIND = RKIND) :: mean_iteration, mean_iteration_loc
    
CONTAINS
    

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
        !La viscosité est calculée pour correspondre au nombre de reynolds choisi
        viscosity = u_c*l_c/re
        
        !Scheme choisi pour termes convectifs
        READ(10, *)
        READ(10, *) scheme

        READ(10, *)
        READ(10, *)
        READ(10, *)
        !Permet de générer des solides à l'intérieur de l'écoulement
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
    
    
    
    !subroutine pour l'ecriture des données dans un fichier
    SUBROUTINE write_output_file(iteration)
    IMPLICIT NONE
        
        INTEGER(KIND = IKind), INTENT(IN) :: iteration
        
        CHARACTER(LEN = StrLen) :: name
        INTEGER(KIND = RKind) :: i, j
        
        !$OMP BARRIER
        
        !$OMP SINGLE
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
        !$OMP END SINGLE
        
    END SUBROUTINE write_output_file
    
    
    !Permet de rapidement tester la valeur de certaines variables
    SUBROUTINE debug(iteration)
    IMPLICIT NONE
        
        INTEGER(KIND = IKind), INTENT(IN) :: iteration
        
        CHARACTER(LEN = StrLen) :: name
        INTEGER(KIND = RKind) :: i, j
        
        !$OMP BARRIER
        
        !$OMP SINGLE
        WRITE(name, '(I0)') iteration
        name = 'debug/a_opti_' // TRIM(name) // '.dat'
        
        !Ouverture du fichier a ecrire
        OPEN(11, FILE = name)
        
        DO i = 1, n_x*n_y
            WRITE(11, '(A, I3, A, 5(F7.1))') 'i = ', i, ' | a_opti(i) = ', a_opti(i, :)
        END DO
        
        !Fermeture du fichier
        CLOSE(11)
        !$OMP END SINGLE
        
    END SUBROUTINE Debug
    
        
    
    
    !maj à l'Etape 7, 2D
    !subroutine de creation du maillage spatial, contient les coordonnées exactes de chaque pt, en 2D
    SUBROUTINE create_space_grid()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j, k
        INTEGER, DIMENSION(:, :), ALLOCATABLE :: borders_grid
        REAL(KIND = RKIND) :: x, y, poly1, poly2, squares1, squares2
        
        !Allocation des tableaux liés au maillage spatial
        ALLOCATE(space_grid%x(n_x))
        ALLOCATE(space_grid%y(n_y))
        !Les tableaux suivant permettent de gérer des géométries internes plus complexes (cercles, carres, polynomes)
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
        
        !Génère les bordures
        borders_grid(:, :) = 0
        borders_grid(-1:0, :) = 1
        borders_grid(n_x+1:n_x+2, :) = 1
        borders_grid(:, -1:0) = 1
        borders_grid(:, n_y+1:n_y+2) = 1
        
        !Attribue les gradients aux surfaces des bords (calcul du vecteur normal à la surface pour pouvoir appliquer les BC de pression)
        space_grid%grad_x(:, :) = 0.0_RKind
        space_grid%grad_y(:, :) = 0.0_RKind
        
        space_grid%grad_y(:, 0) = 1.0_RKind
        space_grid%grad_y(:, n_y+1) = -1.0_RKind
        space_grid%grad_x(0, :) = 1.0_RKind
        space_grid%grad_x(n_x+1, :) = -1.0_RKind
        
        !Gère la création des solides à l'intérieur de l'écoulement en fonction des conditions choisies (fonctions)
        !Calcul des normales aux surfaces
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
        
        !Tableau définissant la positions des bords proches
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
        
        !j\i     -2  -1  0   +1  +2
            !+2  . | . |128| . | . 
            !+1  . | . | 8 | . | . 
            ! 0  16| 1 | x | 2 | 32
            !-1  . | . | 4 | . | . 
            !-2  . | . | 64| . | . 
        
        DEALLOCATE(borders_grid)
        
    END SUBROUTINE create_space_grid
    
    
    
    !Initialisation aux conditions initiales, elles s'appliquent aux deux composantes u et v et à la pression
    SUBROUTINE init_solution()
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j
        
        ALLOCATE(u(n_x, n_y))
        ALLOCATE(v(n_x, n_y))
        ALLOCATE(p(n_x, n_y))
        
        !Assignation pour chaque position en fonction des conditions choisies
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
        
        !Si le point est dans un solide, vitesse et pression nulle (la pression en ces points n'est pas utilisée dans la résolution des systèmes)
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
    
    
    
    !Création et remplissage de la matrice A et allocation des tableaux nécessaires à Jacobi
    SUBROUTINE init_a()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j, k, k_max
        REAL(KIND = RKIND) :: inv_x_2, inv_y_2

        !calcul de la dimension du vecteur solution
        k_max = n_x*n_y
        
        !Allocation de tous les tableaux nécessaires au calcul du champ de pression
        ALLOCATE(a_opti(k_max, 5))
        ALLOCATE(p_vec(k_max))
        ALLOCATE(p_vec_temp(k_max))
        ALLOCATE(residual(k_max))
        ALLOCATE(b(k_max))
        ALLOCATE(conjugate(k_max))
        ALLOCATE(a_mul_conj(k_max))
        ALLOCATE(a_mul_residual(k_max))
        
        
        
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
        
        
        !La matrice A ne contiendrait que 5 diagonale, par soucis d'optimisation on la reduit donc au tableau ci-dessous (5 colonnes, et k_max lignes)
        !Pour une grille 200*200 la matrice A pèserait 13Go sur la RAM en double précision, il est donc important de libérer cet espace
        !De plus ce choix permet de vectoriser facilement les calculs matriciels
        
        a_opti(:, :) = 0.0_RKind
        
        !L'implementation ci-dessous fonctionne pour des surfaces non plane, c'est pourquoi l'on utilise les vecteurs normaux (ou gradients) aux surfaces
        DO j = 1, n_y
            DO i = 1, n_x

                k = (j - 1)*n_x + i
                
                IF (space_grid%borders(i,j) < 0) THEN
                
                    a_opti(k, 3) = 0_RKind
                    
                !Conditions limites de dérivée nulle
                ELSE IF (MOD(space_grid%borders(i,j), 2) == 1) THEN
                    
                    a_opti(k, 3) = -ABS(space_grid%grad_x(i-1, j))/dx - ABS(space_grid%grad_y(i-1, j))/dy
                    a_opti(k, 4) = ABS(space_grid%grad_x(i-1, j))/dx
                    
                    IF (ABS(space_grid%grad_y(i-1, j)) < 1E-14_RKind) THEN
                        CONTINUE
                    ELSE IF (space_grid%grad_y(i-1, j) > 0) THEN
                        a_opti(k, 5) = ABS(space_grid%grad_y(i-1, j))/dy
                    ELSE
                        a_opti(k, 1) = ABS(space_grid%grad_y(i-1, j))/dy
                    END IF

                ELSE IF (MOD(space_grid%borders(i,j), 4)/2 == 1) THEN
                    
                    a_opti(k, 3) = -ABS(space_grid%grad_x(i+1, j))/dx - ABS(space_grid%grad_y(i+1, j))/dy
                    a_opti(k, 2) = ABS(space_grid%grad_x(i+1, j))/dx
                    
                    IF (ABS(space_grid%grad_y(i+1, j)) < 1E-14_RKind) THEN
                        CONTINUE
                    ELSE IF (space_grid%grad_y(i+1, j) > 0) THEN
                        a_opti(k, 5) = ABS(space_grid%grad_y(i+1, j))/dy
                    ELSE
                        a_opti(k, 1) = ABS(space_grid%grad_y(i+1, j))/dy
                    END IF
                    
                ELSE IF (MOD(space_grid%borders(i,j), 8)/4 == 1) THEN
                    
                    a_opti(k, 3) = -ABS(space_grid%grad_x(i, j-1))/dx - ABS(space_grid%grad_y(i, j-1))/dy
                    a_opti(k, 5) = ABS(space_grid%grad_y(i, j-1))/dy
                    
                    IF (ABS(space_grid%grad_x(i, j-1)) < 1E-14_RKind) THEN
                        CONTINUE
                    ELSE IF (space_grid%grad_x(i, j-1) > 0) THEN
                        a_opti(k, 4) = ABS(space_grid%grad_x(i, j-1))/dx
                    ELSE
                        a_opti(k, 2) = ABS(space_grid%grad_x(i, j-1))/dx
                    END IF
                    
                ELSE IF (MOD(space_grid%borders(i,j), 16)/8 == 1) THEN
                    
                    a_opti(k, 3) = -ABS(space_grid%grad_x(i, j+1))/dx - ABS(space_grid%grad_y(i, j+1))/dy
                    a_opti(k, 1) = ABS(space_grid%grad_y(i, j+1))/dy
                    
                    IF (ABS(space_grid%grad_x(i, j+1)) < 1E-14_RKind) THEN
                        CONTINUE
                    ELSE IF (space_grid%grad_x(i, j+1) > 0) THEN
                        a_opti(k, 4) = ABS(space_grid%grad_x(i, j+1))/dx
                    ELSE
                        a_opti(k, 2) = ABS(space_grid%grad_x(i, j+1))/dx
                    END IF
                    
                ELSE 
                    
                    !sinon on applique les coefficients de l'équation
                    a_opti(k, 3) = -2_RKind*(inv_x_2 + inv_y_2)
                    a_opti(k, 4) = inv_x_2
                    a_opti(k, 2) = inv_x_2
                    a_opti(k, 5) = inv_y_2
                    a_opti(k, 1) = inv_y_2
                END IF
            END DO
        END DO
        
    END SUBROUTINE init_a
    
    
    
    !maj à l'étape 3, 2D
    !Réalise toute les initialisations
    SUBROUTINE initialisation()
    
    IMPLICIT NONE
        
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
        
        !Matrice pour le calcul de la pression
        CALL init_a()
        
        i = 0
        !Ecriture de la solution initiale
        CALL write_output_file(i)
        
    END SUBROUTINE initialisation
    
    
    
    !calcul du pas de temps pour une itération, cfl imposé constant, mais dépend de u ou de v 
    SUBROUTINE compute_time_step()

    IMPLICIT NONE

        INTEGER(KIND = IKIND) :: i, j
        REAL(KIND = RKIND) :: dt_min = 0, dt_temp = 0

        !On cherche le dt le plus contraignant, càd le plus petit

        !$OMP BARRIER
        
        !$OMP SINGLE
        !initialisation du min en utilisant la condition du fourrier
        dt_min = fo/viscosity*dx**2

        !recherche des u_max, v_max
        DO i = 1, n_x
            DO j = 1, n_y
                
                !calcul du minimum potentiel sur v
                dt_temp = cfl*dx/ABS(u(i,j))

                IF (dt_temp < dt_min) THEN
                    dt_min = dt_temp
                END IF
                
                !calcul du minimum potentiel sur v
                dt_temp = cfl*dy/ABS(v(i,j))

                IF (dt_temp < dt_min) THEN
                    dt_min = dt_temp
                END IF
            END DO
        END DO
        
        dt = dt_min
        !$OMP END SINGLE
        
    END SUBROUTINE compute_time_step
    
    
    
    !remplit le vecteur b (terme de droite de l'equation de poisson)
    SUBROUTINE fill_b ()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: k_max = 0, i, j, k
        
        !$OMP SINGLE
        k_max = n_x*n_y
        !$OMP END SINGLE
        
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j, k)
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
                    
                END IF
                
            END DO
        END DO
        !$OMP END DO
        
    END SUBROUTINE fill_b
    
    
    
    !Calcul de la norme d'un vecteur
    SUBROUTINE norm_2(vec, norm)
    
    IMPLICIT NONE
        
        REAL(KIND = Rkind), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: vec
        
        INTEGER(KIND = IKIND) :: vec_size = 0, i
        REAL(KIND = RKIND) :: norm
        
        !$OMP SINGLE
        vec_size = SIZE(vec, 1)
        !$OMP END SINGLE
        
        ! $OMP WORKSHARE
        !$OMP SINGLE
        norm = SUM(vec(:)**2_RKind)
        !$OMP END SINGLE
        ! $OMP END WORKSHARE
        
        !$OMP SINGLE
        norm = SQRT(norm)
        !$OMP END SINGLE
        
    END SUBROUTINE norm_2
    
    
    !Integrale de la pression, et correction pour la maintenir nulle
    SUBROUTINE pressure_integral_correction(integral)
    
    IMPLICIT NONE
        
        REAL(KIND = RKind) :: integral
        INTEGER(KIND = IKIND) :: i, j
        
        !$OMP BARRIER
        integral = 0_RKIND
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j) REDUCTION(+:integral)
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
        !$OMP END DO
        !$OMP SINGLE
        integral = integral/(l_x*l_y)
        p_vec(:) = p_vec(:) - integral
        !$OMP END SINGLE
        
    END SUBROUTINE pressure_integral_correction
    
    
    
    ! VOIR AVEC LE PROF SI TEMPS, probleme d'allocation
    ! FUNCTION a_opti_mul(vec)
    
    ! IMPLICIT NONE
        
    !     REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: a_opti_mul
        
        
    ! END FUNCTION a_mul
    
    
    
    !Methode de Jacobi (resolution iterative de systeme lineaire)
    SUBROUTINE jacobi_method()
    
    IMPLICIT NONE
        
        !point d'arret de jacobi
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001
        !Arret forcé de jacobi
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000
        REAL(KIND = RKind) :: upper_norm = 1, lower_norm = 1, integral = 0
        INTEGER(KIND = RKind) :: i, j, k_max = 0, iteration = 0
        
        !$OMP BARRIER
        
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        !$OMP END DO
        
        !$OMP SINGLE
        k_max = n_x*n_y
        
        !Calcul du résidu sous forme vectorisée
        residual(:) = - b(:)
        
        residual(1:k_max) = residual(1:k_max) + a_opti(1:k_max, 3)*p_vec(1:k_max)
        residual(2:k_max) = residual(2:k_max) + a_opti(2:k_max, 2)*p_vec(1:k_max-1)
        residual(1:k_max-1) = residual(1:k_max-1) + a_opti(1:k_max-1, 4)*p_vec(2:k_max)
        residual(n_x+1:k_max) = residual(n_x+1:k_max) + a_opti(n_x+1:k_max, 1)*p_vec(1:k_max-n_x)
        residual(1:k_max-n_x) = residual(1:k_max-n_x) + a_opti(1:k_max-n_x, 5)*p_vec(1+n_x:k_max)
        
        lower_norm = 1
        upper_norm = 1
        
        iteration = 0
        !$OMP END SINGLE
        
        CALL pressure_integral_correction(integral)
        
        !$OMP BARRIER
        
        DO WHILE (upper_norm/lower_norm > RTol)
            
            !$OMP SINGLE
            !Amélioration de la solution
            p_vec_temp(:) = p_vec(:)
            !$OMP END SINGLE
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
            DO i = 1, k_max
                IF (ABS(a_opti(i, 3)) > 1E-15_RKind) THEN
                    p_vec(i) = p_vec_temp(i) - residual(i)/a_opti(i, 3)
                END IF
            END DO
            !$OMP END DO
            
            CALL pressure_integral_correction(integral)
            
            !$OMP BARRIER
            
            CALL norm_2(p_vec, lower_norm)
            
            !$OMP BARRIER
            
            !$OMP SINGLE
            p_vec_temp(:) = p_vec(:)-p_vec_temp(:)
            !$OMP END SINGLE
            
            CALL norm_2(p_vec_temp, upper_norm)
            
            !$OMP BARRIER
            
            !$OMP SINGLE
            !Calcul du résidu
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
            !$OMP END SINGLE
            
        END DO
        
        !$OMP BARRIER
        
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
        !Récupération de la pression dans le tableau 2D
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        !$OMP END DO
        
        !$OMP SINGLE
        PRINT*, 'Jacobi :', iteration, ' iterations | integrale(p) = ', integral
        
        !Pour le benchmark
        mean_iteration_loc = iteration
        !$OMP END SINGLE
        
    END SUBROUTINE jacobi_method
    
    
    
    !Methode de Gauss-Siedel (resolution iterative de systeme lineaire)
    SUBROUTINE gauss_siedel_method()
    
    IMPLICIT NONE
        
        !point d'arret de la méthode
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001
        !Arret forcé de la méthode
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        
        k_max = n_x*n_y
        
        iteration = 0
        
        CALL pressure_integral_correction(integral)
        
        lower_norm = 1
        upper_norm = 1
        DO WHILE (upper_norm/lower_norm > RTol)
            
            !Amélioration de la solution
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
            p_vec_temp(:) = p_vec(:)-p_vec_temp(:)
            CALL norm_2(p_vec_temp, upper_norm)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max de Gauss-Siedel atteint (', IterationMax, ' )'
                STOP
            END IF
            
        END DO
        
        !Récupération de la solution dans le tableau 2D
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        PRINT*, 'Gauss-Siedel :', iteration, ' iterations | integrale(p) = ', integral
        
        !Pour le benchmark
        mean_iteration_loc = iteration
        
    END SUBROUTINE gauss_siedel_method
    
    
    
    !Methode de Surrelaxation (resolution iterative de systeme lineaire)
    SUBROUTINE successive_over_relaxation_method()
    
    IMPLICIT NONE
        
        !point d'arret de la méthode
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001
        !Arret forcé de la méthode
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        REAL(KIND = RKind) :: sor_coeff
        
        !Coefficient de surrelaxation (successive over-relaxation)
        !Choix arbitraire dépendant du problème traité (0<coeff<2)
        sor_coeff = 1.3
        
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        
        k_max = n_x*n_y
        
        iteration = 0
        
        CALL pressure_integral_correction(integral)
        
        lower_norm = 1
        upper_norm = 1
        DO WHILE (upper_norm/lower_norm > RTol)
            
            !Amélioration de la solution
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
            p_vec_temp(:) = p_vec(:)-p_vec_temp(:)
            CALL norm_2(p_vec_temp, upper_norm)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max de Surrelaxation atteint (', IterationMax, ' )'
                STOP
            END IF
            
        END DO
        
        !Récupération de la solution sous forme 2D
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        PRINT*, 'Surrelaxation :', iteration, ' iterations | integrale(p) = ', integral
        
        !Benchmark
        mean_iteration_loc = iteration
        
    END SUBROUTINE successive_over_relaxation_method
    
    
    
    !Methode de la descente du gradient (resolution iterative de systeme lineaire)
    SUBROUTINE steepest_gradient_method()
    
    IMPLICIT NONE
        
        !point d'arret de la méthode
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001
        !Arret forcé de la méthode
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000
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
        
        k_max = n_x*n_y
        
        
        iteration = 0
        
        CALL pressure_integral_correction(integral)
        
        !Calcul du résidu
        residual(:) = b(:)
        residual(1:k_max) = residual(1:k_max) - a_opti(1:k_max, 3)*p_vec(1:k_max)
        residual(2:k_max) = residual(2:k_max) - a_opti(2:k_max, 2)*p_vec(1:k_max-1)
        residual(1:k_max-1) = residual(1:k_max-1) - a_opti(1:k_max-1, 4)*p_vec(2:k_max)
        residual(n_x+1:k_max) = residual(n_x+1:k_max) - a_opti(n_x+1:k_max, 1)*p_vec(1:k_max-n_x)
        residual(1:k_max-n_x) = residual(1:k_max-n_x) - a_opti(1:k_max-n_x, 5)*p_vec(1+n_x:k_max)
        
        !Choix de la première direction
        conjugate(:) = residual(:)
        
        lower_norm = 1
        upper_norm = 1
        DO WHILE (upper_norm/lower_norm > RTol)
            
            p_vec_temp(:) = p_vec(:)
            
            !A(:,:)*conjugué(:)
            a_mul_conj(:) = 0.0_RKind
            a_mul_conj(1:k_max) = a_mul_conj(1:k_max) + a_opti(1:k_max, 3)*conjugate(1:k_max)
            a_mul_conj(2:k_max) = a_mul_conj(2:k_max) + a_opti(2:k_max, 2)*conjugate(1:k_max-1)
            a_mul_conj(1:k_max-1) = a_mul_conj(1:k_max-1) + a_opti(1:k_max-1, 4)*conjugate(2:k_max)
            a_mul_conj(n_x+1:k_max) = a_mul_conj(n_x+1:k_max) + a_opti(n_x+1:k_max, 1)*conjugate(1:k_max-n_x)
            a_mul_conj(1:k_max-n_x) = a_mul_conj(1:k_max-n_x) + a_opti(1:k_max-n_x, 5)*conjugate(1+n_x:k_max)
            
            !Choix de la distance à parcourir dans la direction de conjugate
            alpha = DOT_PRODUCT(conjugate, residual)/DOT_PRODUCT(conjugate, a_mul_conj)
            !Mise à jour du vecteur solution
            p_vec(:) = p_vec_temp(:) + alpha*conjugate(:)
            
            !Calcul du résidu
            residual(:) = b(:)
            residual(1:k_max) = residual(1:k_max) - a_opti(1:k_max, 3)*p_vec(1:k_max)
            residual(2:k_max) = residual(2:k_max) - a_opti(2:k_max, 2)*p_vec(1:k_max-1)
            residual(1:k_max-1) = residual(1:k_max-1) - a_opti(1:k_max-1, 4)*p_vec(2:k_max)
            residual(n_x+1:k_max) = residual(n_x+1:k_max) - a_opti(n_x+1:k_max, 1)*p_vec(1:k_max-n_x)
            residual(1:k_max-n_x) = residual(1:k_max-n_x) - a_opti(1:k_max-n_x, 5)*p_vec(1+n_x:k_max)
            
            !Choix de la nouvelle direction
            conjugate(:) = residual(:)
            
            CALL pressure_integral_correction(integral)
            
            CALL norm_2(p_vec, lower_norm)
            p_vec_temp(:) = p_vec(:)-p_vec_temp(:)
            CALL norm_2(p_vec_temp, upper_norm)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max du Gradient conjugue atteint (', IterationMax, ' )'
                STOP
            END IF
            
        END DO
        
        !Récupération du vecteur solution
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        PRINT*, 'Descente du gradient :', iteration, ' iterations | integrale(p) = ', integral
        
        !Benchmark
        mean_iteration_loc = iteration
        
    END SUBROUTINE steepest_gradient_method
    
    
    
    !Methode du gradient conjugué (méthode itérative)
    SUBROUTINE conjugate_gradient_method()
    
    IMPLICIT NONE
        
        !point d'arret de la méthode
        REAL(KIND = RKind), PARAMETER :: RTol = 0.001
        !Arret forcé de la méthode
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 20000
        REAL(KIND = RKind) :: upper_norm, lower_norm, time1, time2, convergence, integral, r_norm0, r_norm
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        REAL(KIND = RKind) :: alpha, beta, r_r
        
        !Tentative initiale
        DO j = 1, n_y
            DO i = 1, n_x
                p_vec((j-1)*n_x + i) = p(i, j)
            END DO
        END DO
        
        k_max = n_x*n_y
        
        
        iteration = 0
        
        !CALL pressure_integral_correction(integral)
        
        !Calcul du résidu
        residual(:) = b(:)
        residual(1:k_max) = residual(1:k_max) - a_opti(1:k_max, 3)*p_vec(1:k_max)
        residual(2:k_max) = residual(2:k_max) - a_opti(2:k_max, 2)*p_vec(1:k_max-1)
        residual(1:k_max-1) = residual(1:k_max-1) - a_opti(1:k_max-1, 4)*p_vec(2:k_max)
        residual(n_x+1:k_max) = residual(n_x+1:k_max) - a_opti(n_x+1:k_max, 1)*p_vec(1:k_max-n_x)
        residual(1:k_max-n_x) = residual(1:k_max-n_x) - a_opti(1:k_max-n_x, 5)*p_vec(1+n_x:k_max)
        
        !Choix de la direction
        conjugate(:) = residual(:)
        
        lower_norm = 1
        upper_norm = 1
        
        !CALL norm_2(residual, upper_norm)
        !lower_norm = upper_norm
        
        DO WHILE (upper_norm/lower_norm > RTol)
            
            p_vec_temp(:) = p_vec(:)
            
            !A(:,:)*conjugué(:)
            a_mul_conj(:) = 0.0_RKind
            a_mul_conj(1:k_max) = a_mul_conj(1:k_max) + a_opti(1:k_max, 3)*conjugate(1:k_max)
            a_mul_conj(2:k_max) = a_mul_conj(2:k_max) + a_opti(2:k_max, 2)*conjugate(1:k_max-1)
            a_mul_conj(1:k_max-1) = a_mul_conj(1:k_max-1) + a_opti(1:k_max-1, 4)*conjugate(2:k_max)
            a_mul_conj(n_x+1:k_max) = a_mul_conj(n_x+1:k_max) + a_opti(n_x+1:k_max, 1)*conjugate(1:k_max-n_x)
            a_mul_conj(1:k_max-n_x) = a_mul_conj(1:k_max-n_x) + a_opti(1:k_max-n_x, 5)*conjugate(1+n_x:k_max)
            
            
            !Choix de la distance, la simplification classique de ce calcul ne fonctionne pas ici à cause des modifications liées à l'intégrale nulle
            alpha = DOT_PRODUCT(conjugate, residual)/DOT_PRODUCT(conjugate, a_mul_conj)
            !Mise à jour du vecteur solution
            p_vec(:) = p_vec_temp(:) + alpha*conjugate(:)
            
            !Calcul du résidu
            residual(:) = b(:)
            residual(1:k_max) = residual(1:k_max) - a_opti(1:k_max, 3)*p_vec(1:k_max)
            residual(2:k_max) = residual(2:k_max) - a_opti(2:k_max, 2)*p_vec(1:k_max-1)
            residual(1:k_max-1) = residual(1:k_max-1) - a_opti(1:k_max-1, 4)*p_vec(2:k_max)
            residual(n_x+1:k_max) = residual(n_x+1:k_max) - a_opti(n_x+1:k_max, 1)*p_vec(1:k_max-n_x)
            residual(1:k_max-n_x) = residual(1:k_max-n_x) - a_opti(1:k_max-n_x, 5)*p_vec(1+n_x:k_max)
            
            
            a_mul_residual(:) = 0.0_RKind
            a_mul_residual(1:k_max) = a_mul_residual(1:k_max) + a_opti(1:k_max, 3)*residual(1:k_max)
            a_mul_residual(2:k_max) = a_mul_residual(2:k_max) + a_opti(2:k_max, 2)*residual(1:k_max-1)
            a_mul_residual(1:k_max-1) = a_mul_residual(1:k_max-1) + a_opti(1:k_max-1, 4)*residual(2:k_max)
            a_mul_residual(n_x+1:k_max) = a_mul_residual(n_x+1:k_max) + a_opti(n_x+1:k_max, 1)*residual(1:k_max-n_x)
            a_mul_residual(1:k_max-n_x) = a_mul_residual(1:k_max-n_x) + a_opti(1:k_max-n_x, 5)*residual(1+n_x:k_max)
            
            !Choix du décalage vis-à-vis du résidu
            !La simplification classique de ce calcul ne fonctionne pas ici à cause des modifications liées à l'intégrale nulle
            beta = DOT_PRODUCT(residual, a_mul_conj)/DOT_PRODUCT(conjugate, a_mul_conj)
            !beta = -DOT_PRODUCT(conjugate, p_vec_temp(:))/DOT_PRODUCT(conjugate, a_mul_conj)
            !Mise à jour de la direction
            conjugate(:) = residual(:) + beta*conjugate(:)
            
            !Rectification via l'intégrale
            CALL pressure_integral_correction(integral)
            
            !CALL norm_2(residual, upper_norm)
            !PRINT*, upper_norm
            
            CALL norm_2(p_vec, lower_norm)
            a_mul_residual(:) = p_vec(:)-p_vec_temp(:)
            CALL norm_2(a_mul_residual, upper_norm)
            
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max du Gradient conjugue atteint (', IterationMax, ' )'
                STOP
            END IF
            
        END DO
        
        !Récupération de la pression sous la forme d'un tableau 2D
        DO j = 1, n_y
            DO i = 1, n_x
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        PRINT*, 'Gradient conjugue :', iteration, ' iterations | integrale(p) = ', integral
        
        !Benchmark
        mean_iteration_loc = iteration
        
    END SUBROUTINE conjugate_gradient_method
    
    
    !Schéma centré d'ordre 2 pour l'estimation de la vitesse
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
    
    
    
    !Estimation du profil de vitesse pour une itération temporelle
    SUBROUTINE speed_guess()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i, j
        !permet de déterminer la direction du upwind (1 si backward, 0 sinon)
        INTEGER(KIND = IKIND) :: upwind_x, upwind_y
        
        !$OMP BARRIER
        
        ! $OMP WORKSHARE
        !$OMP SINGLE
        !Conditions limites
        u_temp(1, :) = setup%u(1)
        u_temp(n_x, :) = setup%u(2)
        u_temp(:, 1) = setup%u(3)
        u_temp(:, n_y) = setup%u(4)
        
        v_temp(1, :) = setup%v(1)
        v_temp(n_x, :) = setup%v(2)
        v_temp(:, 1) = setup%v(3)
        v_temp(:, n_y) = setup%v(4)
        !$OMP END SINGLE
        ! $OMP END WORKSHARE
        
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
        DO j = 1, n_y
            DO i = 1, n_x
                IF (space_grid%borders(i, j) < 0) THEN
                    u(i, j) = 0_RKind
                    v(i, j) = 0_RKind
                    p(i, j) = 0_RKind
                END IF
            END DO
        END DO
        !$OMP END DO
        
        
        
        !Calcul avec scheme régressif upwind d'ordre 1
        IF (scheme == 'UR1') THEN
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j, upwind_x, upwind_y)
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
            !$OMP END DO
        
        !Calcul avec scheme centré d'ordre 4
        ELSE IF (scheme == 'CD4') THEN
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i)
            DO i = 2, n_x-1
                CALL cd_scheme_2(i, 2)
                CALL cd_scheme_2(i, n_y - 1)
            END DO
            !$OMP END DO
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(j)
            DO j = 3, n_y-2
                CALL cd_scheme_2(2, j)
                CALL cd_scheme_2(n_x-1, j)
            END DO
            !$OMP END DO
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(j)
            DO j = 3, n_y-2
                    
                !Calcul de u
                u_temp(3:n_x-2,j) = u(3:n_x-2,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind)) &
                + u(2:n_x-3,j)*dt*(viscosity/dx**2_RKind + 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + u(3:n_x-2,j-1)*dt*(viscosity/dy**2_RKind + 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                + u(4:n_x-1,j)*dt*(viscosity/dx**2_RKind - 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + u(3:n_x-2,j+1)*dt*(viscosity/dy**2_RKind - 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                - u(1:n_x-4,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                - u(3:n_x-2,j-2)*dt*v(3:n_x-2,j)/(12_RKind*dy) &
                + u(5:n_x-1,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                + u(3:n_x-2,j+2)*dt*v(3:n_x-2,j)/(12_RKind*dy)
                
                !Calcul de v
                v_temp(3:n_x-2,j) = v(3:n_x-2,j)*(1.0 - 2.0*viscosity*dt*(1.0_RKind/dx**2_RKind + 1.0_RKind/dy**2_RKind)) &
                + v(2:n_x-3,j)*dt*(viscosity/dx**2_RKind + 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + v(3:n_x-2,j-1)*dt*(viscosity/dy**2_RKind + 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                + v(4:n_x-1,j)*dt*(viscosity/dx**2_RKind - 4_RKind*u(3:n_x-2,j)/(6_RKind*dx)) &
                + v(3:n_x-2,j+1)*dt*(viscosity/dy**2_RKind - 4_RKind*v(3:n_x-2,j)/(6_RKind*dy)) &
                - v(1:n_x-4,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                - v(3:n_x-2,j-2)*dt*v(3:n_x-2,j)/(12_RKind*dy) &
                + v(5:n_x-1,j)*dt*u(3:n_x-2,j)/(12_RKind*dx) &
                + v(3:n_x-2,j+2)*dt*v(3:n_x-2,j)/(12_RKind*dy)
                
            END DO
            !$OMP END DO
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
            !Schéma centré d'ordre 2 pour les bords
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    IF (space_grid%borders(i, j)/16 >= 0) THEN
                        CALL cd_scheme_2(i, j)
                    END IF
                END DO
            END DO
            !$OMP END DO
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    IF (MOD(space_grid%borders(i, j), 16) /= 0) THEN
                        u_temp(i, j) = 0_RKIND
                        v_temp(i, j) = 0_RKIND
                    END IF
                END DO
            END DO
            !$OMP END DO
            
        !Calcul avec scheme centré d'ordre 2
        ELSE IF (scheme == 'CD2') THEN
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    CALL cd_scheme_2(i, j)
                END DO
            END DO
            
            !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
            DO j = 2, n_y-1
                DO i = 2, n_x-1
                    IF (MOD(space_grid%borders(i, j), 16) /= 0) THEN
                        u_temp(i, j) = 0_RKIND
                        v_temp(i, j) = 0_RKIND
                    END IF
                END DO
            END DO
            !$OMP END DO
            
        ELSE
            
            !$OMP SINGLE
            PRINT*, 'Mauvaise définition du schema (fichier input.dat)'
            STOP
            !$OMP END SINGLE
            
        END IF
        
        !$OMP BARRIER
        
    END SUBROUTINE speed_guess


    !Résout l'équation de poisson contenant les prédictions de vitesse, pour une itération temporelle
    SUBROUTINE compute_pressure()
    
    IMPLICIT NONE
        
        !$OMP BARRIER
        
        !remplissage de b à chaque itération temporelle, là où A est constant
        CALL fill_b()
        
        !$OMP BARRIER
        
        !resolution de l'equation de poisson avec une méthode itérative
        
        CALL jacobi_method()
        !CALL gauss_siedel_method()
        !CALL successive_over_relaxation_method()
        !CALL steepest_gradient_method()
        !CALL conjugate_gradient_method()
        
        !$OMP BARRIER
        
    END SUBROUTINE compute_pressure
    
    
    
    !Corrige la vitesse de façon à satisfaire l'équation de continuité grâce au calcul de la pression
    SUBROUTINE adjust_speed()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, j
        
        !$OMP BARRIER
        
        ! $OMP WORKSHARE
        !$OMP SINGLE
        !Conditions limites
        u(1, :) = setup%u(1)
        u(n_x, :) = setup%u(2)
        u(:, 1) = setup%u(3)
        u(:, n_y) = setup%u(4)
        
        v(1, :) = setup%v(1)
        v(n_x, :) = setup%v(2)
        v(:, 1) = setup%v(3)
        v(:, n_y) = setup%v(4)
        !$OMP END SINGLE
        ! $OMP END WORKSHARE
        
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
        DO j = 1, n_y
            DO i = 1, n_x
                IF (space_grid%borders(i, j) < 0) THEN
                    u(i, j) = 0_RKind
                    v(i, j) = 0_RKind
                    p(i, j) = 0_RKind
                END IF
            END DO
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(j)
        DO j = 2, n_y-1
            u(2:n_x-1, j) = u_temp(2:n_x-1, j) - ((p(3:n_x, j) - p(1:n_x-2, j))/(2.0_RKind*dx))*dt/density
            v(2:n_x-1, j) = v_temp(2:n_x-1, j) - ((p(2:n_x-1, j+1) - p(2:n_x-1, j-1))/(2.0_RKind*dy))*dt/density
        END DO
        !$OMP END DO
        
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i, j)
        DO j = 2, n_y-1
            DO i = 2, n_x-1
                IF (MOD(space_grid%borders(i, j), 16) /= 0) THEN
                    u(i, j) = 0_RKIND
                    v(i, j) = 0_RKIND
                END IF
            END DO
        END DO
        !$OMP END DO
        
    END SUBROUTINE adjust_speed
    
    
    
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
        
        CALL omp_set_num_threads(1)
        !$OMP PARALLEL DEFAULT(SHARED)
        
        !Boucle temporelle du calcul
        DO WHILE ((last_iteration .EQV. .FALSE.) .AND. (i < NMax))
            
            !$OMP BARRIER
            
            CALL compute_time_step()
            
            !$OMP BARRIER
            
            IF (t + dt > t_f) THEN
                last_iteration = .TRUE.
                EXIT
            END IF
            
            !$OMP BARRIER
            
            !$OMP SINGLE
            PRINT*, 'dt = ', dt
            !$OMP END SINGLE
            
            !$OMP BARRIER
            
            !appelle le calcul pour cette itération
            CALL speed_guess()
            
            !$OMP BARRIER
            
            CALL compute_pressure()
            
            !$OMP BARRIER
            
            !$OMP SINGLE
            IF (i == 0) THEN
                mean_iteration = mean_iteration + mean_iteration_loc
            END IF
            !$OMP END SINGLE
            
            CALL adjust_speed()
            
            !$OMP BARRIER
            
            !$OMP SINGLE
            i = i + 1
            t = t + dt
            !$OMP END SINGLE
            
            !écrit dans un fichier toute les frames
            IF (MOD(i, frame) == 0) THEN
                CALL write_output_file(i)
            END IF
            
            !$OMP SINGLE
            PRINT*, 't = ', t
            PRINT*, '-------------------'
            !$OMP END SINGLE
            
        END DO
        
        !$OMP END PARALLEL
        
        DEALLOCATE(u_temp)
        DEALLOCATE(v_temp)
        
    END SUBROUTINE resolution_loop
    
END MODULE global



PROGRAM main

USE global

IMPLICIT NONE
    
    REAL(KIND = RKIND) :: time1, time2, mean_time
    INTEGER, DIMENSION(5) :: mesh_size
    INTEGER :: i, j, nb_tests
    CHARACTER(LEN = StrLen) :: name
    
    !Taille des maillages pour le benchmark
    mesh_size(1) = 21
    mesh_size(2) = 31
    mesh_size(3) = 51
    mesh_size(4) = 101
    mesh_size(5) = 201
    
    !nb de calculs à chaque maillage
    nb_tests = 3
    
    
    ! !Programme de benchmark
    ! OPEN(12, FILE = 'benchmark/relaxation_short_5.dat')
    
    ! DO i = 1, SIZE(mesh_size)
        
    !     mean_time = 0.0_RKind
    !     mean_iteration = 0.0_RKIND
        
    !     DO j = 1, nb_tests
        
    !         CALL CPU_TIME(time1)
            
    !         !récupération des données du problème
    !         name = 'input.dat'
    !         CALL read_input_file(name)
            
    !         n_x = mesh_size(i)
    !         n_y = mesh_size(i)
    !         t_f = 0.5
            
    !         CALL initialisation()
            
    !         CALL resolution_loop()
            
            
            
    !         DEALLOCATE(space_grid%x)
    !         DEALLOCATE(space_grid%y)
    !         DEALLOCATE(space_grid%borders)
    !         DEALLOCATE(u)
    !         DEALLOCATE(v)
    !         DEALLOCATE(a_opti)
    !         DEALLOCATE(b)
    !         DEALLOCATE(p)
    !         DEALLOCATE(p_vec)
    !         DEALLOCATE(p_vec_temp)
    !         DEALLOCATE(residual)
    !         DEALLOCATE(space_grid%grad_x)
    !         DEALLOCATE(space_grid%grad_y)
    !         DEALLOCATE(conjugate)
    !         DEALLOCATE(a_mul_conj)
    !         DEALLOCATE(a_mul_residual)
            
    !         CALL CPU_TIME(time2)
            
    !         PRINT*, 'mesh = ', n_x, ' | j = ', j
    !         PRINT*, '-------------------------'
            
    !         mean_time = time2 - time1
            
    !     END DO
        
    !     mean_time = mean_time/REAL(nb_tests, RKind)
    !     mean_iteration = mean_iteration/REAL(nb_tests, RKind)
        
    !     WRITE(12, *) dx, mean_iteration, mean_time
    ! END DO
    
    ! CLOSE(12)
    
    
    !Programme classique
    CALL CPU_TIME(time1)
    
    !récupération des données du problème
    name = 'input.dat'
    CALL read_input_file(name)
    
    CALL initialisation()
    
    CALL resolution_loop()
    
    
    !Désallocation des variables
    DEALLOCATE(space_grid%x)
    DEALLOCATE(space_grid%y)
    DEALLOCATE(space_grid%borders)
    DEALLOCATE(u)
    DEALLOCATE(v)
    DEALLOCATE(a_opti)
    DEALLOCATE(b)
    DEALLOCATE(p)
    DEALLOCATE(p_vec)
    DEALLOCATE(p_vec_temp)
    DEALLOCATE(residual)
    DEALLOCATE(space_grid%grad_x)
    DEALLOCATE(space_grid%grad_y)
    DEALLOCATE(conjugate)
    DEALLOCATE(a_mul_conj)
    DEALLOCATE(a_mul_residual)
    
    CALL CPU_TIME(time2)
    
    PRINT*, 'temps de resolution = ', time2 - time1
    
END PROGRAM main