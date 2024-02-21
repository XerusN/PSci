MODULE global

IMPLICIT NONE

    !definition des tailles des reels(RKind) et des entiers(IKind)
    INTEGER, PARAMETER :: RKind = SELECTED_REAL_KIND(10,200)
    INTEGER, PARAMETER :: IKind = SELECTED_INT_KIND(5)

    !Constantes de dimension des tableaux
    INTEGER, PARAMETER :: Nx = 201
    
    !Constante de taille pour le nom des fichiers
    INTEGER, PARAMETER :: StrLen = 40
    
    
    
    !-----------------------------------------------------
    
    
    
    !Variables d'entree

    !nombre d'itérations temporelles
    INTEGER(KIND = IKind) :: n_t
    
    !dimension spatiale du maillage [m]
    REAL(KIND = Rkind) :: l

    !position initiale du centre du pic de la gaussienne [m]
    REAL(KIND = Rkind) :: x_0 

    !vitesse de convection [m/s]
    REAL(KIND = Rkind) :: c

    !temps total d'observation du phénomène [s]
    REAL(KIND = Rkind) :: t_f

    !CFL nombre de courant adimensionnel
    REAL(KIND = Rkind) :: cfl

    !écart type [m]
    REAL(KIND = Rkind) :: delta

    !conditions limites a gauche et a droite [m/s] (même unité que u)
    REAL(KIND = Rkind) :: boundary_condition_left, boundary_condition_right
    
    !definition de la vitesse [m/s] comme var globale, tableau de scalaire pour la partie 1, on ne stocke pas les itérations
    REAL(KIND = RKind), DIMENSION(Nx) :: u
    
    !definition du maillage spatial, exprime les coordonnées en [m] de la case
    REAL(KIND = RKind), DIMENSION(Nx) :: space_grid

    !definition du pas spatial [m] et du pas de temps [s]
    REAL(KIND = RKind) :: dx, dt

    
    
CONTAINS
    


    !maj à l'Etape 1, 1D
    SUBROUTINE read_input_file(name)
    IMPLICIT NONE
        
        CHARACTER(LEN = StrLen), INTENT(IN) :: name
        
        !Ouverture du fichier
        OPEN(10, FILE = name)
        
        
        !Saut des deux premieres lignes
        READ(10, *)
        READ(10, *)
        
        !longueur espace
        READ(10, *) l
        
        !Nb d'iteration en temps
        READ(10, *) n_t
        !Temps de fin
        READ(10, *) t_f
        
        !Vitesse initiale
        READ(10, *) c
        !CFL
        READ(10, *) cfl
        
        !Variables d'initialisations
        READ(10, *) x_0
        READ(10, *) delta
        
        !Conditions aux limites spatiales
        READ(10, *) boundary_condition_left
        READ(10, *) boundary_condition_right
        
        
        !Fermeture du fichier
        CLOSE(10)
        
    END SUBROUTINE read_input_file
    
    
    
    !maj à l'Etape 1, 1D
    !subroutine pour l'ecriture des données dans un fichier
    SUBROUTINE write_output_file(name)
    IMPLICIT NONE
        
        CHARACTER(LEN = StrLen), INTENT(IN) :: name
        
        INTEGER(KIND = RKind) :: i
        
        !Ouverture du fichier a ecrire
        OPEN(11, FILE = name)
        
        
        !Ecriture pour chaque position
        DO i = 1, Nx
            WRITE(11, *) space_grid(i), u(i)
        END DO
        
        
        !Fermeture du fichier
        CLOSE(11)
        
    END SUBROUTINE write_output_file
        
    
    
    !maj à l'Etape 1, 1D
    !subroutine de creation du maillage spatial, contient les coordonnées exactes de chaque pt, en 1D
    SUBROUTINE create_space_grid()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i

        !calcul du pas spatial
        dx = l / (Nx - 1)

        !assignation des coordonnées
        DO i = 1, Nx
            space_grid(i) = dx * (i-1)
        END DO

    END SUBROUTINE create_space_grid
    
    
    
    !maj à l'Etape 1, 1D
    !Initialisation aux conditions initiales
    SUBROUTINE init_solution()
    IMPLICIT NONE
        
        INTEGER(KIND = IKind) :: i
        
        !Calcul pour chaque position
        DO i= 1, Nx
            u(i) = EXP(- ((space_grid(i) - x_0)/delta)**2)
        END DO
        
    END SUBROUTINE init_solution
    
    
    
    !maj à l'Etape 1, 1D
    !calcul du profil de vitesse pour une itération temporelle
    SUBROUTINE convection_lineaire_1D()

    IMPLICIT NONE

        INTEGER(KIND = IKind) :: i !compteur
        REAL(KIND = RKind), DIMENSION(Nx) :: u_temp !profil de vitesse stocké temporairement

        !definition des conditions limites a gauche et a droite, inchangées
        u(1) = boundary_condition_left
        u(Nx) = boundary_condition_right
        
        !calcul du n+1
        DO i = 2, Nx-1
            u_temp(i) = u(i) * (1 - c * dt / dx) + u(i-1) * (c * dt /dx)
        END DO

        !assignation du n+1
        DO i = 2, Nx-1
            u(i) = u_temp(i)
        END DO

    END SUBROUTINE convection_lineaire_1D

END MODULE global






PROGRAM main

USE global

IMPLICIT NONE

    INTEGER(KIND = IKind) :: i

    CHARACTER(LEN = StrLen) :: str
    
    !récupération des données du problème
    str = 'input.dat'
    CALL read_input_file(str)
    
    !création du maillage spatial
    CALL create_space_grid()
    
    !calcul du pas de temps à partir du CFL
    dt = cfl * dx / c
    
    !initalisation de la solution grâce aux conditions initiales
    CALL init_solution()
    
    !Ecriture de la solution initiale
    str = 'output/initial_solution.dat'
    CALL write_output_file(str)
    
    
    !Boucle temporelle du calcul
    DO i = 1, n_t
        CALL convection_lineaire_1D()
    END DO
    
    
    !Ecriture de la solutions à t_f
    str = 'output/final_solution.dat'
    CALL write_output_file(str)

END PROGRAM main