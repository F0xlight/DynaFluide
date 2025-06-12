#---------------------------------------------------------------
# Simulation de mécanique des fluides v0.1
#---------------------------------------------------------------
# On supposera le problème invariant selon un axe de sorte à 
# le réduire à deux dimensions et le fluide un gaz parfait 
# incompressible sans viscosité
#---------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

duree = 5       # Durée de l'expérience en s
delta_t = 0.1   # Pas d'intégration en s
delta_l = 1     # Pas de la grille en mm
rho = 1000      # Masse volumique du fluide en g/m^3
p0 = 10**5      # Pression initiale en Pa
R = 8.3         # Constante du gaz parfait en J/mol/K
M = 0.029       # Masse molaire du gaz en kg/mol
g = 9.81        # Acceleration de la pesanteur en m/s^2

taillex = 100   # Taille de la fenetre de simulation en mm
tailley = 100

Nx, Ny = taillex//delta_l, tailley//delta_l
Nt = int(duree/delta_t)

# [ux,uy,p,modifiable] ; (0,0) en haut à gauche
grille = np.array([[[[0,0,p0,True]for _ in range(Ny)] for _ in range(Nx)] for _ in range(Nt)], dtype = float)

# Conditions initiales 

for t in range(Nt) :
    for x in range(Nx) :
        grille[t,x,0,0] = 10 
        grille[t,x,0,3] = False


# Schéma de résolution

for t in range(1,Nt) :
    for x in range(1,Nx) :
        for y in range(1,Ny) :
            epsilon = 1e-8  # Small value to avoid division by zero
            grille[t,x,y,0] = grille[t-1,x,y,0] - delta_l*((grille[t-1,x,y,2]-grille[t-1,x-1,y,2])/(rho*delta_l))/(grille[t-1,x,y,0] + epsilon) - grille[t-1,x,y,1]*(grille[t-1,x,y,0]-grille[t-1,x,y-1,0])/(grille[t-1,x,y,0] + epsilon)
            grille[t,x,y,1] = grille[t-1,x,y,1] + delta_l*(rho*delta_l**3*g - (grille[t-1,x,y,2]-grille[t-1,x,y-1,2])/(rho*delta_l))/(grille[t-1,x,y,1] + epsilon) - grille[t-1,x,y,0]*(grille[t-1,x,y,1]-grille[t-1,x-1,y,1])/(grille[t-1,x,y,1] + epsilon)
            
    # Affichage
    plt.clf()
    data = grille[t,:,:,0]  # Use ux as colormap
    plt.imshow(data.T, origin='lower', cmap='viridis', extent=[0, Nx, 0, Ny])
    plt.colorbar(label='ux')
    plt.title(f"t = {t*delta_t:.2f} s")
    plt.pause(0.01)

