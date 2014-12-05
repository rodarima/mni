# encoding: utf-8

import numpy as np
import pyqtgraph as pg

from math import sin, pi

# -------------------------------------------------------------------

def f01 (x):

    return np.sin (np.pi + x)

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

#n = int(input('\n\n	Introduce el número de nodos:   '))
n=5

x = np.linspace(xini, xfin, n)
h = x[1] - x[0]
h2i = 1.0 / (h*h)

# -------------------------------------------------------------------
# Matriz de rigidez.

a = np.zeros((n,n))
for i in range(1,n-1):			# Excluimos el primero y último nodos.
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

#plt.figure ()

tt = np.linspace (xini, xfin, 101)
yy = analitica (tt)
#plt.plot (tt,yy, 'r-', lw=2)
#pg.plot(tt, yy)

win = pg.GraphicsWindow(title="d1")
#win.resize(1000,600)
#win.setWindowTitle('pyqtgraph example: Plotting')

# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)

#p1 = win.addPlot(title="Solución", y=u, x=x, pen=(200,200,200), 
#symbolBrush=(255,0,0), symbolPen='w')
#p1.showGrid(x=True, y=True)
#
#p2 = win.addPlot(title="Analítica", y=yy, x=tt)
#p2.showGrid(x=True, y=True)
#p2.setXLink(p1)
#p2.setYLink(p1)

p1 = win.addPlot(title="Solución", y=u, x=x, pen=(200,200,200), symbolPen='w')
p1.plot(x=tt, y=yy, pen=(200,0,0))


if __name__ == '__main__':
    import sys
    if sys.flags.interactive != 1 or not hasattr(QtCore, 'PYQT_VERSION'):
        pg.QtGui.QApplication.exec_()

#plt.plot (x,u, 'bo')
#plt.grid (True)

# -------------------------------------------------------------------

#plt.show ()

# -------------------------------------------------------------------

