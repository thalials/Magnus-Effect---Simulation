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
import pylab as p
import mpl_toolkits.mplot3d.axes3d as p3
 	
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
            Cd = 0.25 
            
            D = -ro*(v**2)*A*Cd/2
            Dx = D*vx/v
            Dy = D*vy/v
            Dz = D*vz/v
            
            M = -ro*(v**2)*A*Cl/2
            Mx = M*(vz-vy)/v
            My = M*vx/v
            Mz = -M*vx/v
            
            Rx = Dx + Mx
            Ry = Dy + My - m*g
            Rz = Dz + Mz
            
            dvxdt = Rx/m
            dvydt = Ry/m
            dvzdt = Rz/m


        elif v == 0:
            D = 0
            M = 0
            Ry = -m*g
            dvxdt = 0
            dvydt = Ry/m
            dvzdt = 0

        
    return dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt

def eq_dif2(lista_equacoes, tempo):
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
            Cl = 0 # o coeficiente n depende de w pq esse valor já tem o w incluso (valor experimental)
            Cd = 0.25 
            
            D = -ro*(v**2)*A*Cd/2
            Dx = D*vx/v
            Dy = D*vy/v
            Dz = 0
            
            M = -ro*(v**2)*A*Cl/2
            Mx = M*(vz-vy)/v
            My = M*vx/v
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
lista_tempo = np.arange(0, 15, dt)

# Condicoes iniciais: x, y, vx, vy respectivamente
v = 76 # velocidade de lançamento (m/s)
theta = 17 # angulo de lançamento (graus)
condicoes_iniciais = [0, 0, 0, v*cos(theta*pi/180), v*sin(theta*pi/180), 0]


# Solucoes
solucao2 = odeint(eq_dif2, condicoes_iniciais, lista_tempo) # somente c om Arrasto
condicoes_iniciais = [0, 0, 0, v*cos(theta*pi/180), v*sin(theta*pi/180),0]
solucao = odeint(eq_dif, condicoes_iniciais, lista_tempo) # com efeito Magnus e Arrasto
x = solucao[:,0]
y = solucao[:,1]
z = solucao[:,2]
vx = solucao[:,2]
vy = solucao[:,3]
vz = solucao[:,4]

x1 = solucao2[:,0]
y1 = solucao2[:,1]
vx1 = solucao2[:,2]
vy1 = solucao2[:,3]

##Graficos:

plt.plot(x, y, label = '{:.0f}'.format(theta)) 
plt.plot(x1, y1, label = 'Sem Efeito Magnus') #nesse caso devemos usar as equações da fisica normal 
plt.title("Trajetória")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
#plt.axis([0,250,0,70])
plt.grid(True)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


fig=p.figure()
ax = p3.Axes3D(fig)
ax.plot(solucao[:,0], solucao[:,2], solucao[:,1], label = "Com Efeito Magnus")
ax.plot(solucao2[:,0], solucao2[:,2], solucao2[:,1],'r', label = "Sem Efeito Magnus")
ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")
ax.set_zlabel("y (m)")
ax.set_title("Trajetória da bola")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))