#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""

from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt

# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Nelectrons = 200  # change this to have more or fewer atoms
Nions = 16
dt = 1E-5  # pas d'incrémentation temporel

simulation_steps = 1000

# Déclaration de variables physiques "Typical values"
DIM = 2 #Nombre de degrés de liberté de la simulation 
mass = 4E-3/6E23 # helium mass
Relectron = 0.005 # just for visualisation purposes
Rion = 0.1
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature

#### CANEVAS DE FOND ####
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
animation = canvas(width=750, height=500)
animation.range = L

#### ARÊTES DE BOÎTE 2D ####
d = L/2+Relectron
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d,-d,0), vector(d,-d,0), vector(d,d,0), vector(-d,d,0), vector(-d,-d,0)])

#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
electrons = [] # Objet qui contiendra les sphères pour l'animation
p = [] # quantité de mouvement des sphères
electrons_pos = [] # position des sphères
pavg = sqrt(2*mass*(DIM/2)*k*T) # average kinetic energy in 3D p**2/(2mass) = (3/2)kT : Principe de l'équipartition de l'énergie en thermodynamique statistique classique

for i in range(Nelectrons):
    x = L*random()-L/2 # position aléatoire qui tient compte que l'origine est au centre de la boîte
    y = L*random()-L/2
    z = 0
    electrons.append(simple_sphere(pos=vector(x,y,z), radius=Relectron, color=gray))
    electrons_pos.append(vec(x,y,z)) # liste de la position initiale de toutes les sphères
#    theta = pi*random() # direction de coordonnées sphériques, superflue en 2D
    phi = 2*pi*random() # direction aléatoire pour la quantité de mouvement
    px = pavg*cos(phi)  # qte de mvt initiale selon l'équipartition
    py = pavg*sin(phi)
    pz = 0
    p.append(vector(px,py,pz)) # liste de la quantité de mouvement initiale de toutes les sphères

ions = []
ions_pos = []

for i in range(Nions):
    grid_size = int(sqrt(Nions))
    x = (i % grid_size) * L / grid_size - L / 2 + L / (2 * grid_size)
    y = (i // grid_size) * L / grid_size - L / 2 + L / (2 * grid_size)
    z = 0
    ions.append(simple_sphere(pos=vector(x,y,z), radius=Rion, color=color.red))
    ions_pos.append(vec(x,y,z))

#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions():
    hitlist = []   # initialisation
    for i in range(Nelectrons):
        ai = electrons_pos[i]
        for j in range(Nions) :
            aj = ions_pos[j]
            dr = ai - aj   # la boucle dans une boucle itère pour calculer cette distance vectorielle dr entre chaque paire de sphère
            if dr.mag <= Rion:   # test de collision où mag2(dr) qui retourne la norme élevée au carré de la distance intersphère dr
                hitlist.append([i,j]) # liste numérotant toutes les paires de sphères en collision
    return hitlist

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre

particule_a_suivre = 0
particule_a_suivre_parcours = []

for step in range(simulation_steps):
    rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!

    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    vitesse = []   # vitesse instantanée de chaque sphère
    deltax = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for i in range(Nelectrons):
        vitesse.append(p[i]/mass)   # par définition de la quantité de nouvement pour chaque sphère
        deltax.append(vitesse[i] * dt)   # différence avant pour calculer l'incrément de position
        electrons[i].pos = electrons_pos[i] = electrons_pos[i] + deltax[i]  # nouvelle position de l'atome après l'incrément de temps dt

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES MURS DE LA BOÎTE ####
    for i in range(Nelectrons):
        loc = electrons_pos[i]
        if abs(loc.x) > L/2:
            if loc.x < 0: p[i].x =  abs(p[i].x)  # renverse composante x au mur de gauche
            else: p[i].x =  -abs(p[i].x)   # renverse composante x au mur de droite
        if abs(loc.y) > L/2:
            if loc.y < 0: p[i].y = abs(p[i].y)  # renverse composante y au mur du bas
            else: p[i].y =  -abs(p[i].y)  # renverse composante y au mur du haut

    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist = checkCollisions()

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for ij in hitlist:

        # définition de nouvelles variables pour chaque paire de sphères en collision
        i = ij[0]  # extraction du numéro des 2 sphères impliquées à cette itération
        j = ij[1]
        #ptot = p[i]+p[j]   # quantité de mouvement totale des 2 sphères
        #mtot = 2*mass    # masse totale des 2 sphères
        #Vcom = ptot/mtot   # vitesse du référentiel barycentrique/center-of-momentum (com) frame
        posi = electrons_pos[i]   # position de chacune des 2 sphères
        posj = ions_pos[j]
        vi = p[i]/mass   # vitesse de chacune des 2 sphères
        #vj = p[j]/mass
        rrel = posi-posj  # vecteur pour la distance entre les centres des 2 sphères
        #vrel = vj-vi   # vecteur pour la différence de vitesse entre les 2 sphères

        # exclusion de cas où il n'y a pas de changements à faire
        if vi.mag2 == 0: continue  # exactly same velocities si et seulement si le vecteur vrel devient nul, la trajectoire des 2 sphères continue alors côte à côte
        #if rrel.mag > Relectron: continue  # one atom went all the way through another, la collision a été "manquée" à l'intérieur du pas deltax

        #if particule_a_suivre in ij:
        #    particule_a_suivre_parcours.append((
        #        electrons_pos[particule_a_suivre],
        #        (step+1)*dt
        #    ))

        perp = vector(-rrel.y, rrel.x) # perpendiculaire à la collision
        random_direction = -pi*random()
        new_direction = cos(random_direction) * perp + sin(random_direction) * rrel

        p[i] = p[i].mag * new_direction

        # calcule la distance et temps d'interpénétration des sphères dures qui ne doit pas se produire dans ce modèle
        #dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta) où theta is the angle between vrel and rrel:
        #dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
        #alpha = asin(dy/(2*Relectron))  # alpha is the angle of the triangle composed of rrel, path of atom j, and a line from the center of atom i to the center of atom j where atome j hits atom i
        #d = (2*Relectron)*cos(alpha)-dx # distance traveled into the atom from first contact
        #deltat = d/vrel.mag         # time spent moving from first contact to position inside atom

        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        #posi = posi-vi*deltat   # back up to contact configuration
        #posj = posj-vj*deltat
        #pcomi = p[i]-mass*Vcom  # transform momenta to center-of-momentum (com) frame
        #pcomj = p[j]-mass*Vcom
        #rrel = hat(rrel)    # vecteur unitaire aligné avec rrel
        #pcomi = pcomi-2*dot(pcomi,rrel)*rrel # bounce in center-of-momentum (com) frame
        #pcomj = pcomj-2*dot(pcomj,rrel)*rrel
        #p[i] = pcomi+mass*Vcom # transform momenta back to lab frame
        #p[j] = pcomj+mass*Vcom
        #electrons_pos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time, ramenant au même temps où sont rendues les autres sphères dans l'itération
        #electrons_pos[j] = posj+(p[j]/mass)*deltat

