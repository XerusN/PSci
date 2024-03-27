MODULE global
!Poisson
IMPLICIT NONE

    !definition des tailles des reels(RKind) et des entiers(IKind)
    INTEGER, PARAMETER :: RKind = SELECTED_REAL_KIND(10,200)
    INTEGER, PARAMETER :: IKind = SELECTED_INT_KIND(10)
    
    !Constante de taille pour le nom des fichiers
    INTEGER, PARAMETER :: StrLen = 40
    
    
    
    !-----------------------------------------------------
    
    
    
    !Variables d'entree
    
    !Taille du maillage spatial
    INTEGER(KIND = IKind) :: n_x, n_y
    
    !dimension spatiale du maillage [m]
    REAL(KIND = RKind) :: l_x, l_y
    
    !stocke la valeur de la pression dans l'espace
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: p
    !matrice A de l'equation de poisson
    REAL(KIND = RKind), DIMENSION(:, :), ALLOCATABLE :: a
    !vecteur b de l'équation de poisson
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: b
    !Stockage de p sous forme vectorielle
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: p_vec, p_vec_temp, jacobi_r
    !stockage temporaire de a
    REAL(KIND = RKind), DIMENSION(:,:), ALLOCATABLE :: a_loc
    
    !structure qui contient les coordonnées en x et y d'un point du maillage
    TYPE COORDS
        REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: x, y
    END TYPE

    !definition du maillage spatial, exprime les coordonnées en [m] de la case
    TYPE(COORDS) :: space_grid
    
    !variable de fréquence d'affichage
    INTEGER(KIND = IKind) :: frame
    
    !definition du pas spatial [m] et du pas de temps [s]
    REAL(KIND = RKind) :: dx, dy

    
    
    
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
        
        !longueur espace
        READ(10, *)
        READ(10, *) l_x, l_y
        
        READ(10, *)
        
        !fréquence d'affichage
        READ(10, *)
        READ(10, *) frame
        
        
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
        WRITE(11, *) 'VARIABLES = "X", "Y", "P"'
        WRITE(11, '(1X, A, ES20.13, A, I4, A, I4, A)') 'ZONE T="', REAL(iteration), &
            '   seconds", I=', n_x, ', J=', n_y, ', DATAPACKING=POINT'
        
        !Ecriture pour chaque position
        DO i = 1, n_x
            DO j = 1, n_y
                WRITE(11, '(3(ES20.13, 1X))') space_grid%x(i), space_grid%y(j), p(i, j)
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
    
    
    
    !maj à l'Etape 4, 2D
    !Initialisation aux conditions initiales, elles s'appliquent aux deux composantes u et v
    SUBROUTINE init_solution()
    IMPLICIT NONE

        
        ALLOCATE(p(n_x, n_y))
        
        !Assignation de base pour chaque position
        p(:, :) = 0

        !Application des conditions limites, dès maintenant pour tracer la solution init
        !CL basse y = 0
        p(:, 1) = 1 - 2*SINH(space_grid%x(:))
        
        !CL haute y = 1
        p(:, n_y) = EXP(space_grid%x(:)) - 2*SINH(space_grid%x(:))

        !CL gauche x = 0
        p(1, :) = 1

        !CL droite x = 1
        p(n_x, :) = EXP(space_grid%y(:)) - 2*SINH(space_grid%x(n_x))
        
    END SUBROUTINE init_solution
    
    
    
    !maj à l'étape 4, 2D
    !Création et remplissage de la matrice A et du vecteur b du système
    SUBROUTINE matrix_fill()
    
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


        !Remplissage de la matrice A
        !initialisation
        a(:, :) = 0

        !calcul des coefficients
        inv_x_2 = REAL(1)/(dx**2.0)
        inv_y_2 = REAL(1)/(dy**2.0)

        
        DO j = 1, n_y
            DO i = 1, n_x

                k = (j - 1)*n_x + i
                
                !les lignes où on a une condition limite ont juste 1 sur la diagonale
                IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN

                    a(k, k) = 1
                    
                    
                ELSE 
                    !sinon on applique les coefficients de l'équation
                    a(k, k) = REAL(-2*(inv_x_2 + inv_y_2))
                    a(k, k + 1) = inv_x_2
                    a(k, k - 1) = inv_x_2
                    a(k, k + n_x) = inv_y_2
                    a(k, k - n_x) = inv_y_2

                END IF

            END DO
        END DO
        
        
    END SUBROUTINE matrix_fill
    
    
    !remplit le vecteur b
    SUBROUTINE fill_b ()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: k_max, i, j, k
        
        k_max = n_x*n_y
        ALLOCATE(b(k_max))
        
        
        DO j = 1, n_y
            DO i = 1, n_x

                k = (j - 1)*n_x + i
                IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN

                    !application de la CL à b
                    !CL basse y = 0
                    IF (j == 1) THEN
                        b(k) = 1 - 2*SINH(space_grid%x(i))
                    !CL haute y = 1
                    ELSE IF (j == n_y) THEN
                        b(k) = EXP(space_grid%x(i)) - 2*SINH(space_grid%x(i))
                    !CL gauche x = 0
                    ELSE IF (i == 1) THEN
                        b(k) = 1
                    !CL droite x = 1
                    ELSE
                        b(k) = EXP(space_grid%y(j)) - 2*SINH(space_grid%x(n_x))
                    END IF
                    
                    
                ELSE 
                    
                    !valeur de b pour tous les autres points
                    b(k) = (space_grid%x(i)**2 + space_grid%y(j)**2)*EXP(space_grid%x(i)*space_grid%y(j))
                    b(k) = b(k) - 2*SINH(space_grid%x(i))

                END IF
            END DO
        END DO
        
        
        
        
    END SUBROUTINE fill_b
    
    
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
        
        !récupération des données du problème
        name = 'input.dat'
        CALL read_input_file(name)
        
        !création du maillage spatial
        CALL create_space_grid()
    
        !initalisation de la solution grâce aux conditions initiales
        CALL init_solution()
        
        

        !Remplissage des matrices de l'equation
        CALL matrix_fill()
        CALL fill_b()
        
        i = 0
        !Ecriture de la solution initiale
        CALL write_output_file(i)
        
        ! i = -3
        ! !Ecriture de la solution initiale
        ! CALL write_output_file(i)
        
    END SUBROUTINE initialisation
    
    
    
    !maj à l'Etape 4, 2D
    !pivot de gauss
    SUBROUTINE gauss_elimination(k_max)

    IMPLICIT NONE
    

        INTEGER(KIND = IKind) :: k_max
        
        INTEGER(KIND = IKind) :: i, j, k !compteur
        REAL(KIND = RKind) :: pivot, pivot_loc, test, time1, time2
        
        CALL CPU_TIME(time1)
        
        a_loc(:,:) = a(:,:)
        
        !Triangularisation
        DO i = 1, k_max-1
            pivot = a_loc(i, i)
            DO j = i+1, k_max
                pivot_loc = a_loc(j, i)/pivot
                a_loc(j, i:) = a_loc(j, i:) - pivot_loc*a_loc(i, i:)
                b(j) = b(j) - pivot_loc*b(i)
            END DO
        END DO
        
        
        !Resolution backward
        p_vec(k_max) = b(k_max)/a_loc(k_max, k_max)
        DO i = k_max-1, 1, -1
            p_vec(i) = (b(i) - SUM(a_loc(i, i+1:)*p_vec(i+1:)))/a_loc(i, i)
        END DO
        
        
        DO i = 1, n_x
            DO j = 1, n_y
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        
        i = -1
        CALL write_output_file(i)
        
        
        CALL CPU_TIME(time2)
        
        PRINT*, 'Gauss for a grid size of ', n_x, ' : ', time2 - time1, ' seconds'

    END SUBROUTINE gauss_elimination
    
    
    
    !Calcul la norme 2 d'un vecteur
    SUBROUTINE norm_2(vec, norm)
    
    IMPLICIT NONE
        
        REAL(KIND = Rkind), DIMENSION(:), ALLOCATABLE :: vec
        
        INTEGER(KIND = IKIND) :: vec_size, i
        REAL(KIND = RKIND) :: norm
        
        vec_size = SIZE(vec, 1)
        
        norm = SUM(vec(:)**2)
        
        
        norm = SQRT(norm)
        
    END SUBROUTINE norm_2
    
    
    !Methode de Jacobi (resolution iterative de systeme lineaire)
    SUBROUTINE jacobi_method()
    
    IMPLICIT NONE
        
        REAL(KIND = RKind), PARAMETER :: RTol = 0.00001     !point d'arret de jacobi
        INTEGER(KIND = IKind), PARAMETER :: IterationMax = 100000       !Arret forcé de jacobi
        REAL(KIND = RKind) :: initial_norm, r_norm, time1, time2
        INTEGER(KIND = RKind) :: i, j, k_max, iteration
        
        CALL CPU_TIME(time1)
        
        !Tentative initiale
        p_vec(:) = 0
        p_vec_temp(:) = 0
        
        k_max = SIZE(p_vec)
        
        DO i = 1, k_max
            jacobi_r(i) = SUM(a(i, :)*p_vec(:)) - b(i)
        END DO
        
        CALL norm_2(jacobi_r, initial_norm)
        r_norm = initial_norm
        iteration = 0
        
        DO WHILE (r_norm > RTol*initial_norm)
            p_vec_temp(:) = p_vec(:)
            DO i = 1, k_max
                p_vec(i) = p_vec_temp(i) - jacobi_r(i)/a(i, i)
            END DO
            !p_vec(:) = p_vec_temp(:) - jacobi_r(:)/a((/(/j, j/), j=1, k_max/))
            
            
            jacobi_r(:) = - b(:)
            DO j = 1, k_max
                jacobi_r(:) = jacobi_r(:) + a(:, j)*p_vec(j)
            END DO
            
            CALL norm_2(jacobi_r, r_norm)
            iteration = iteration + 1
            IF (iteration >= IterationMax) THEN
                PRINT*, 'Nombre d iteration max de jacobi atteint (', IterationMax, ' )'
                EXIT
            END IF
            
            !permet d'écrire les guess itermediaires
            IF (MOD(iteration, frame) == 0) THEN
                DO i = 1, n_x
                    DO j = 1, n_y
                        p(i, j) = p_vec((j-1)*n_x+i)
                    END DO
                END DO
                
                CALL write_output_file(iteration)
            END IF
            
        END DO
        
        CALL CPU_TIME(time2)
        
        CALL write_output_file(iteration)
        
        
        
        PRINT*, 'Jacobi for a grid size of ', n_x, ' : ', time2 - time1, ' seconds (', iteration, ' iterations)'
        
    END SUBROUTINE jacobi_method
    
    
    !Resolution du probleme
    SUBROUTINE solve()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, k_max, j, k
        
        p(:,:) = 0
        
        k_max = n_x*n_y
        
        
        !resolution avec jacobi
        CALL jacobi_method()
        
        
        !Resolution avec gauss
        CALL gauss_elimination(k_max)
        
        
        
        
        
        
        
    END SUBROUTINE solve
    
    
    !Solution analytique
    SUBROUTINE analytical_solving()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, k_max, j, k
        
        DO i = 1, n_x
            DO j = 1, n_y
                p(i, j) = EXP(space_grid%x(i)*space_grid%y(j)) - 2*SINH(space_grid%x(i))
            END DO
        END DO
        ! i = -2
        ! CALL write_output_file(i)
        
    END SUBROUTINE analytical_solving
    
    

END MODULE global






PROGRAM main

USE global

IMPLICIT NONE
    
    !INTEGER(KIND = IKind) :: i, j
    
    CALL initialisation()
    
    CALL analytical_solving()
    
    CALL solve()
    
    
    
    DEALLOCATE(space_grid%x)
    DEALLOCATE(space_grid%y)
    DEALLOCATE(a)
    DEALLOCATE(a_loc)
    DEALLOCATE(b)
    DEALLOCATE(p)
    DEALLOCATE(p_vec)
    DEALLOCATE(p_vec_temp)
    DEALLOCATE(jacobi_r)
    
        
END PROGRAM main