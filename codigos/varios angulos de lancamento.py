# -*- coding: utf-8 -*-
"""
Autoria:
 Cicero Tiago Carneiro Valentim 
 Thalia Loiola Silva
 Kathleen da Silva Nascimento

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
w = 180                 # velocidade angular (rad/s)
g = 9.8                 # Aceleração gravitacional
m = 0.04                # massa da esfera (kg)
n = 1.8*1e-5            # viscosidade dinamica (Pa.s)
I = (2/3)*m*(r**2)      # Momento de inercia da esfera

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
            Cl = 0.12   # o coeficiente n depende de w pq esse valor já tem o w incluso (valor experimental)
            Cd = 0.25 
            
            D = -ro*(v**2)*A*Cd/2
            Dx = D*vx/v
            Dy = D*vy/v
            Dz = D*vz/v
            
            mod = sqrt(2*( (vz-vy)**2 + vx**2 + vx**2)) # modulo de (w x v)
            M = -ro*(v**2)*A*Cl/2
            Mx = M*(vz-vy)/mod
            My = M*vx/mod
            Mz = -M*vx/mod
            
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
            
            mod = sqrt(2*( (vz-vy)**2 + vx**2 + vx**2)) # modulo de (w x v)
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
lista_tempo = np.arange(0, 15, dt)

# Condicoes iniciais: x, y, vx, vy respectivamente
v =400           # velocidade de lançamento (m/s)
theta = 17          # angulo de lançamento (graus)
unidade = ['°']     # Graus

# Solucoes
condicoes_iniciais = [0, 0, 0, v*cos(theta*pi/180), v*sin(theta*pi/180),0]
solucao = odeint(eq_dif, condicoes_iniciais, lista_tempo)   # com efeito Magnus e Arrasto
solucao2 = odeint(eq_dif2, condicoes_iniciais, lista_tempo) # somente com Arrasto

# Com efeito Magnus e Arrasto
x = solucao[:,0]
y = solucao[:,1]
z = solucao[:,2]
vx = solucao[:,3]
vy = solucao[:,4]
vz = solucao[:,5]

# Somente com Arrasto
x1 = solucao2[:,0]
y1 = solucao2[:,1]
vx1 = solucao2[:,3]
vy1 = solucao2[:,4]

# Gerar graficos para varios angulos de lancamento

# Grafico em 2D
lista_de_angulos = []
lista_de_alcances = []

for i in range(10):
    theta = 4 + 7*i 
    condicoes_iniciais = [0, 0,0, v*cos(theta*pi/180), v*sin(theta*pi/180),0]
    solucao = odeint(eq_dif, condicoes_iniciais, lista_tempo)   # com efeito Magnus e Arrasto
    x = solucao[:,0]
    y = solucao[:,1]
    alcance = sqrt((x[-1]**2)+(z[-1]**2))
    lista_de_alcances.append(alcance)
    lista_de_angulos.append(theta)
    plt.plot(x, y,'--', label = '{:.1f}°'.format(theta))
    
plt.title("Trajetória da bola")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.grid(True)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints=1, frameon=False, labelspacing=1, title='Ângulo de lançamento:')
plt.show()

for i in range(len(lista_de_angulos)):
    plt.plot(lista_de_angulos[i],lista_de_alcances[i],'o', label = "{:.1f}°".format(lista_de_angulos[i]))
plt.title("Ângulo versus Alcance")
plt.xlabel("Ângulo (°)")
plt.ylabel("Alcance (m)")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints=0, frameon=False, labelspacing=1, title='Ângulo de lançamento:')
plt.grid(True)
plt.show()

#
#fig=p.figure()
#ax = p3.Axes3D(fig)
#for i in range(10):
#    theta = 4 + 7*i 
#    condicoes_iniciais = [0, 0,0, v*cos(theta*pi/180), v*sin(theta*pi/180),0]
#    solucao = odeint(eq_dif, condicoes_iniciais, lista_tempo)   # com efeito Magnus e Arrasto
#    x = solucao[:,0]
#    y = solucao[:,1]
#    z = solucao[:,2]
#    ax.plot(x,z,y,'--', label = '{:.1f}°'.format(theta))
    
#ax.set_xlabel("x (m)")
#ax.set_ylabel("z (m)")
#ax.set_zlabel("y (m)")
#ax.set_title("Trajetória da bola")
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
