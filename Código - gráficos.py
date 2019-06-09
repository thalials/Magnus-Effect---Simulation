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

#r = 0.121  	
r = 0.1225		# raio da bola (m)
A = pi*(r**2)          	# Área transversal (m²)
V = 4*pi*(r**3)/3   	# Volume da esfera (m³)
ro = 1.17            	# densidade volumétrica (kg/m³)
#Cd = 0.42           # Coeficiente de arrasto
w = 12          # velocidade angular (rad/s)
g = 10               # Aceleração gravitacional
#m = 0.3 
m = 0.62             # massa da esfera (kg)
n = 1.8*1e-5   #viscosidade dinamica (Pa.s)
I = (2/3)*m*(r**2)

def eq_dif( lista_equacoes, tempo):
    x = lista_equacoes[0]
    y = lista_equacoes[1]
    
    vx = lista_equacoes[2]
    vy = lista_equacoes[3]

    if y <= r:
        y = 0
#        w = 0
        dxdt = 0
        dydt = 0
        dvxdt = 0
        dvydt = 0
#        dwdt = 0 
        
    else:
        v = sqrt((vx**2) + (vy**2))
        dxdt = vx
        dydt = vy
       
        if v > 0 :
            Cl = w*r/v
            Re = 2*ro*r*v/n
            Cd = 24/Re + (2.6*(Re/5))/(1 + ((Re/5)**1.52)) + 0.411*((Re/263000)**(-7.94))/(1+((Re/263000)**(-8))) + (Re**0.8)/461000
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
lista_tempo = np.arange(0, 20, dt)

# Condicoes iniciais: x, y, vx, vy, w respectivamente
condicoes_iniciais = [0, 126.5, 0, 0]

solucao = odeint(eq_dif,condicoes_iniciais,lista_tempo)
# Solucoes 
x = solucao[:,0]
y = solucao[:,1]
vx = solucao[:,2]
vy = solucao[:,3]
#w = solucao[:,4]
#print(w[-1])
print(y[-1])

#Graficos:
# 1) Gráfico da trajetória 
#plt.plot(lista_tempo, w)
#plt.title("w versus tempo")
#plt.xlabel('tempo')
#plt.ylabel('w')
#
#plt.grid()
#plt.show()


# Grafico da velocidade em x versus tempo 
plt.plot(x, y)
plt.title("trajetoria")
plt.xlabel("x")
plt.ylabel("y")
#plt.axis([4, 4.1, -0.00001, 0.00001])
plt.grid()
plt.show()
#
#plt.plot(lista_tempo,solucao[:,3])
#plt.title("vy vs tempo")
#plt.xlabel("tempo")
#plt.ylabel("Velocidade vertical")
#plt.grid()
#plt.show()
#
#plt.plot(solucao[:,2],solucao[:,3], '--')
#plt.xlabel("vx")
#plt.ylabel("vy")
#plt.axis("equal")
#plt.grid()
#plt.show()
#
plt.plot(lista_tempo,solucao[:,1])
plt.title("y vs tempo")
plt.xlabel("tempo")
plt.ylabel("y")
plt.grid()
plt.show()