# -*- coding: utf-8 -*-
"""
Autoria:
 Cicero Tiago Carneiro Valentim 
 Thalia Loiola Silva

"""
from math import sqrt, pi, cos, sin
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

def eq_dif(lista_equacoes, tempo):
    x = lista_equacoes[0]
    y = lista_equacoes[1]
    z = lista_equacoes[2]
    
    vx = lista_equacoes[3]
    vy = lista_equacoes[4]
    vz = lista_equacoes[5]

    v = sqrt((vx**2) + (vy**2) + (vz**2))
    
    if y < 0:
        dxdt = 0
        dydt = 0
        dzdt = 0
        dvxdt = 0
        dvydt = 0
        dvzdt = 0
       
    else:        
        dxdt = vx
        dydt = vy
        dzdt = vz
       
        if v > 0 :
            Cl = 0.12 # o coeficiente n depende de w pq esse valor já tem o w incluso (valor experimental)
            Cd = 0 
            
            D = -ro*(v**2)*A*Cd/2
            Dx = D*vx/v
            Dy = D*vy/v
            Dz = 0
            
            mod = sqrt(2*((vz-vy)**2 + vx**2 + vx**2)) # modulo de (w x v)
            M = -ro*(v**2)*A*Cl/2
            Mx = M*(vz-vy)/mod
            My = M*vx/mod
            Mz = 0
            
            Rx = Dx + Mx
            Ry = Dy + My - m*g
            Rz = Dz + Mz
            
            dvxdt = Rx/m
            dvydt = Ry/m
            dvzdt = 0

        elif v == 0:
            D = 0
            M = 0
            Ry = -m*g
            dvxdt = 0
            dvydt = Ry/m
            dvzdt = 0
        
    return dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt
#lista de tempo analisado
dt = 1e-4
lista_tempo = np.arange(0, 4, dt)

# Condicoes iniciais: x, y, vx, vy respectivamente
v = 76              # velocidade de lançamento (m/s)
theta = 17          # angulo de lançamento (graus)

# Solucoes
condicoes_iniciais = [0, 0, 0, v*cos(theta*pi/180), v*sin(theta*pi/180),0]
solucao = odeint(eq_dif, condicoes_iniciais, lista_tempo)   # com efeito Magnus e Arrasto

# Com efeito Magnus e Arrasto
x = solucao[:,0]
y = solucao[:,1]
z = solucao[:,2]
vx = solucao[:,3]
vy = solucao[:,4]
vz = solucao[:,5]

energia_pot = []
energia_c = []
energia_total = []
for i in range(0, len(lista_tempo)):
    x = solucao[:,0]
    y = solucao[:,1]
    z = solucao[:,2]
    vx = solucao[:,3]
    vy = solucao[:,4]
    vz = solucao[:,5]
    
    v1 = sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2)
    
    pot_grav = m*g*y[i]
    en_cinet = m*(v1**2)/2
    
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