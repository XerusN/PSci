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
    REAL(KIND = RKind), DIMENSION(:), ALLOCATABLE :: p_vec
    !stockage temporaire de a
    REAL(KIND = RKind), DIMENSION(:,:), ALLOCATABLE :: a_loc

    !tableau de stockage intermediaire
    
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
        p(:, 1) = 0   
        
        !CL haute y = 1
        p(:, n_y) = EXP(space_grid%x(:))*SIN(space_grid%y(n_y))

        !CL gauche x = 0
        p(1, :) = SIN(space_grid%y(:))

        !CL droite
        p(n_x, :) = EXP(space_grid%x(n_x))*SIN(space_grid%y(:))
        
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
        ALLOCATE(b(k_max))
        ALLOCATE(p_vec(k_max))

        !Remplissage du vecteur b, on assigne les conditions limites après
        b(:) = 0

        !Remplissage de la matrice A
        !initialisation
        a(:, :) = 0

        !calcul des coefficients
        inv_x_2 = REAL(1)/(dx**2.0)
        inv_y_2 = REAL(1)/(dy**2.0)

        
        ! DO j = 1, n_y
        !     DO i = 1, n_x

        !         k = (j - 1)*n_x + i

        !         !les lignes où on a une condition limite ont juste 1 sur la diagonale
        !         IF ((i == 1) .OR. (i == n_x) .OR. (j == 1) .OR. (j == n_y)) THEN

        !             a(k, k) = 1

        !             !application de la CL à b
        !             b(k) = p(i, j)

        !         ELSE 
        !             !sinon on applique les coefficients de l'équation
        !             a(k, k) = REAL(-2*(inv_x_2 + inv_y_2))
        !             a(k, k + 1) = inv_x_2
        !             a(k, k - 1) = inv_x_2
        !             a(k, k + n_x) = inv_y_2
        !             a(k, k - n_x) = inv_y_2

        !         END IF

                
                
        !     END DO
        ! END DO
        
        
        DO k = 1, k_max
            PRINT*, ((k-MOD(k-1, n_x))/n_x) + 1, MOD(k-1, n_x)+1
            !les lignes où on a une condition limite ont juste 1 sur la diagonale
            IF ((MOD(k-1, n_x)+1 == 1) .OR. (MOD(k-1, n_x)+1 == n_x) .OR. &
                (((k-MOD(k-1, n_x)-1)/n_x) + 1 == 1) .OR. (((k-MOD(k-1, n_x)-1)/n_x) + 1 == n_y)) THEN

                a(k, k) = 1

                !application de la CL à b
                b(k) = p(MOD(k-1, n_x)+1, ((k-MOD(k-1, n_x)-1)/n_x) + 1)

            ELSE 
                !sinon on applique les coefficients de l'équation
                a(k, k) = REAL(-2)*(inv_x_2 + inv_y_2)
                a(k, k + 1) = inv_x_2
                a(k, k - 1) = inv_x_2
                a(k, k + n_x) = inv_y_2
                a(k, k - n_x) = inv_y_2

            END IF

        END DO
        
        
    END SUBROUTINE matrix_fill
    
    
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
        
        i = 0
        !Ecriture de la solution initiale
        CALL write_output_file(i)
        
    END SUBROUTINE initialisation
    
    
    
    !maj à l'Etape 4, 2D
    !pivot de gauss
    SUBROUTINE gauss_elimination(k_max)

    IMPLICIT NONE
    

        INTEGER(KIND = IKind) :: k_max
        
        INTEGER(KIND = IKind) :: i, j, k !compteur
        REAL(KIND = RKind) :: pivot, pivot_loc, test
        
        a_loc(:,:) = a(:,:)
        
        DO i = 1, k_max-1
            pivot = 1.0/a_loc(i, i)
            
            DO j = i+1, k_max
                pivot_loc = a_loc(j, i)*pivot
                DO k = i+1, k_max
                    a_loc(j, k) = a_loc(j, k) - pivot_loc*a_loc(i, k)
                END DO
                b(j) = b(j) - pivot_loc*b(i)
            END DO
        END DO
            
        p_vec(k_max) = b(k_max)/a_loc(k_max, k_max)
        DO i = k_max-1, 1, -1
            pivot = 1.0/a(i,i)
            p_vec(i) = 0
            DO k = i, k_max
                p_vec(i) = p_vec(i) + a_loc(i, k)*p_vec(k)
            END DO
            p_vec(i) = (b(i) - p_vec(i))*pivot
        END DO

    END SUBROUTINE gauss_elimination
    
    
    
    SUBROUTINE solve()
    
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i, k_max, j, k
        
        k_max = n_x*n_y
        CALL gauss_elimination(k_max)
        
        DO i = 1, n_x
            DO j = 1, n_y
                p(i, j) = p_vec((j-1)*n_x+i)
            END DO
        END DO
        
        i = 1
        CALL write_output_file(i)
        
    END SUBROUTINE solve
    
    

END MODULE global






PROGRAM main

USE global

IMPLICIT NONE
    
    INTEGER :: i
    
    CALL initialisation()
    
    CALL solve()
    
    DEALLOCATE(space_grid%x)
    DEALLOCATE(space_grid%y)
    DEALLOCATE(a)
    DEALLOCATE(a_loc)
    DEALLOCATE(b)
    DEALLOCATE(p)
    DEALLOCATE(p_vec)
    
        
END PROGRAM main