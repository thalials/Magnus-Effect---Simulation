# -*- coding: utf-8 -*-
"""
Autoria:
 Cicero Tiago Carneiro Valentim 
 Thalia Loiola Silva

"""
from math import sqrt, pi
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
 	
r = 0.0213		        # raio da bola (m)
A = pi*(r**2)          	# Área transversal (m²)
V = 4*pi*(r**3)/3   	# Volume da esfera (m³)
ro = 1.29            	# densidade volumétrica (kg/m³)
w = 180               # velocidade angular (rad/s)
g = 9.8                # Aceleração gravitacional
m = 0.04                # massa da esfera (kg)
n = 1.8*1e-5            # viscosidade dinamica (Pa.s)
I = (2/3)*m*(r**2)

def eq_dif(lista_equacoes, tempo):
    x = lista_equacoes[0]
    y = lista_equacoes[1]
    
    vx = lista_equacoes[2]
    vy = lista_equacoes[3]

    v = sqrt((vx**2) + (vy**2))
    
    if y < 0:
        dxdt = 0
        dydt = 0
        dvxdt = 0
        dvydt = 0
       
    else:
        
        dxdt = vx
        dydt = vy
       
        if v > 0 :
            Cl = 0.126 #o coeficiente n depende de w pq esse valor já tem o w incluso (valor experimental)
            Re = 2*ro*r*v/n
            Cd = 0 #24/Re + (2.6*(Re/5))/(1 + ((Re/5)**1.52)) + 0.411*((Re/263000)**(-7.94))/(1+((Re/263000)**(-8))) + (Re**0.8)/461000
            cos_alfa = vx/v
            sen_alfa = vy/v
            
            Fd = ro*(v**2)*A*Cd/2
            Fm = ro*(v**2)*A*Cl/2
            
            Rx = -Fd*cos_alfa - Fm*sen_alfa
            dvxdt = Rx/m
            
            Ry = -Fd*sen_alfa + Fm*cos_alfa - m*g
            dvydt = Ry/m

        elif v == 0:
            Fd = 0
            Fm = 0
            Ry = -m*g
            dvxdt = 0
            dvydt = Ry/m

        
    return dxdt, dydt, dvxdt, dvydt

#lista de tempo analisado
dt = 1e-4
lista_tempo = np.arange(0, 16.7, dt)

# Condicoes iniciais: x, y, vx, vy respectivamente
condicoes_iniciais = [0, 0, 130, 50]

solucao = odeint(eq_dif, condicoes_iniciais, lista_tempo)

# Solucoes 
x = solucao[:,0]
y = solucao[:,1]
vx = solucao[:,2]
vy = solucao[:,3]

energia_pot = []
energia_c = []
energia_total = []
for i in range(0, len(lista_tempo)):
    x1 = solucao[:,0][i]
    y1 = solucao[:,1][i]
    vx1 = solucao[:,2][i]
    vy1 = solucao[:,3][i]
    
    v = sqrt(vx1**2 + vy1**2)
    
    pot_grav = m*g*y1
    en_cinet = m*(v**2)/2
    
    energia_c.append(en_cinet)
    energia_pot.append(pot_grav)
    
    energia_total.append(pot_grav + en_cinet)
  
plt.plot(lista_tempo, energia_pot, label ='Energia Potencial Gravitacional')
plt.plot(lista_tempo, energia_c, label = 'Energia Cinética')
plt.plot(lista_tempo, energia_total, label = 'Energia Total')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.title("Conservação da Energia")

plt.xlabel('Tempo (s)')
plt.ylabel('Energia (J)')

plt.show()