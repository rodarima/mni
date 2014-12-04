# encoding: utf-8

from math import sin, pi
import numpy as np
#import pylab as plt
import matplotlib.pyplot as plt

# -------------------------------------------------------------------

def f01 (x):

    return np.sin (pi + x)

# -------------------------------------------------------------------

def f02 (x):

    return 2.0

# -------------------------------------------------------------------

def f03 (x):

    return x * (1.0 - x)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Datos: dominio, función de segundo miembro y solución analítica.

analitica = None

# Caso 1.

'''
xini = 0.0
xfin = pi
second = f01
dir01 = 0.0
dir02 = 0.0
analitica = f01
'''

# Caso 2.

'''
xini = 0.0
xfin = 1.0
second = f02
dir01 = 0.0
dir02 = 0.0
analitica = f03
'''

# Caso 3.

xini = 0.0
xfin = 2.0
second = f02
dir01 = 0.0
dir02 = -2.0
analitica = f03

# -------------------------------------------------------------------
# Mallado.

n = int(input('\n\n    Introduce el número de nodos:   '))

x = np.linspace(xini, xfin, n)
h = x[1] - x[0]
h2i = 1.0 / (h*h)

# -------------------------------------------------------------------
# Matriz de rigidez.

a = np.zeros((n,n))
for i in range(1,n-1):            # Excluimos el primero y último nodos.
    a[i,i] = 2.0
for i in range(1,n-1):
    a[i,i+1] = -1.0
for i in range(2,n-1):
    a[i,i-1] = -1.0

a *= h2i

a[0,0] = 1.0
a[n-1,n-1] = 1.0

# -------------------------------------------------------------------
# Vector segundo miembro.

b = np.zeros (n)
b[0] = dir01
for i in range (1,n-1):
    b[i] = second (x[i])
b[n-1] = dir02

# -------------------------------------------------------------------
# Resolución del sistema de ecuaciones.

u = np.linalg.solve (a,b)

# Añadimos las condiciones de contorno.

#u = list (u)
#u.insert (0, dir01)
#u.append (dir02)

# -------------------------------------------------------------------
# Representación gráfica.

plt.figure ()

if (analitica):
    tt = np.linspace (xini, xfin, 101)
    yy = analitica (tt)
    plt.plot (tt,yy, 'r-', lw=2)

plt.plot (x,u, 'bo')
plt.grid (True)

# -------------------------------------------------------------------

plt.show ()

# -------------------------------------------------------------------

