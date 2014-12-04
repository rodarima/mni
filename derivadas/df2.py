# encoding: utf-8

from math import sin, pi

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np
import matplotlib.pyplot as plt
import time
#import graf01


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Condiciones de contorno.

def bound (x,y,t):

    return 0.0

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Condición inicial.

def init (x,y):

    return x * y * (1.0 - x) * (1.0 - y)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Función de segundo miembro.

def f (x,y,t):

    upx = 1.0 + x
    upy = 1.0 + y
    umx = 1.0 - x
    umy = 1.0 - y
    um2x = 1.0 - 2.0 * x
    um2y = 1.0 - 2.0 * y

#    z = x * y * umx * umy - (t + 1.0) * (um2x * y * umy - 2.0 * y * upx * umy   \
#             + upy * um2x * um2y - x*x * umx * um2y + (4. + x*y) * 2.0 * x * umx)

#    z = x * y * umx * umy + 2.0 * y * umy * (t + 1.0) + 2.0 * x * umx * (t + 1.0)

    z = 0.0

    return z

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Coeficiente de difusión.

def diff (x,y):

    z = np.zeros ((2,2), 'f')

    z[0,0] = 1.
    z[0,1] = 0.2
    z[1,0] = 0.8
    z[1,1] = 0.3

#    z = np.eye (2,2)

#    z [0,0] = 1.0 + x
#    z [0,1] = 1.0 + y
#    z [1,1] = 4.0 + x*y

    return z

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Derivadas del coeficiente de difusión.

def ddiff (x,y):

    ddx = np.zeros ((2,2), 'f')
    ddy = np.zeros ((2,2), 'f')

#    ddx [0,0] = 1.0
#    ddx [1,1] = y

#    ddy [0,1] = 1.0
#    ddy [1,1] = x

    return [ddx,ddy]

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Mallado.

#tfin = int(input ('\n\n    Introduce el instante final:                '))
#nt = int(input ('    Introduce el número de pasos temporales:    '))
#n = int(input ('\n\n    Introduce el número de nodos en abscisas:   '))
#m = int(input ('    Introduce el número de nodos en ordenadas:  '))

tfin = 10
nt = 10
n = 10
m = 10

xini = 0.0
xfin = 1.0
yini = 0.0
yfin = 1.0

nnod = (n + 2) * (m + 2)

x = np.linspace (xini, xfin, n+2)
y = np.linspace (yini, yfin, m+2)
dx = x[1] - x[0]
dy = y[1] - y[0]
dt = tfin / nt

xx,yy=np.meshgrid(x,y)

#xx = np.zeros (nnod, 'f')
#yy = np.zeros (nnod, 'f')
#for j in range (m+2):
#    for i in range (n+2):
#        k = j * (n+2) + i
#        xx[k] = x[i]
#        yy[k] = y[j]

# -------------------------------------------------------------------
# Condición inicial.

uold = np.zeros (nnod, 'f')
for j in range (m+2):
    for i in range (n+2):
        k = j * (n+2) + i
        uold [k] = init (x[i], y[j])
#graf01.graf (xx,yy,uold,0,'Cond. inicial')

# -------------------------------------------------------------------
# Bucle temporal.

plt.ion ()
fig = plt.figure()


for it in range (nt):
    t = (it + 1) * dt

    # Matriz de rigidez.

    a = np.zeros ((nnod,nnod), 'f')
    b = np.zeros (nnod, 'f')

    for j in range (1,m+1):
        for i in range (1,n+1):
            coef = diff (x[i], y[j])
            [dcdx, dcdy] = ddiff (x[i], y[j])

            k = j * (n+2) + i
            no = k + (n+2)
            so = k - (n+2)
            ea = k + 1
            we = k - 1
            ne = no + 1
            nw = no - 1
            se = so + 1
            sw = so - 1
            a[k,k] = 1.0 + (dt/dx) * (dcdx[0,0] + dcdy[1,0]) + (dt/dx) * (dcdx[0,1] + dcdy[1,1])  \
                         + (2.0*dt/dx**2) * coef[0,0] + (2.0*dt/dy**2) * coef[1,1]
            a[k,no] = -(dt/dx) * (dcdx[0,1] + dcdy[1,1]) - (dt/dy**2) * coef[1,1]
            a[k,so] = -(dt/dy**2) * coef[1,1]
            a[k,ea] = -(dt/dx) * (dcdx[0,0] + dcdy[1,0]) - (dt/dx**2) * coef[0,0]
            a[k,we] = -(dt/dx**2) * coef[0,0]

            aux = dt * (coef[0,1] + coef[1,0]) / (4.0 * dx * dy)
            a[k,nw] = aux
            a[k,ne] = aux
            a[k,se] = -aux
            a[k,sw] = -aux

            b[k] = uold [k] + dt * f(x[i], y[j],t)

    # Condiciones de contorno.

    j = 0                         # Frontera inferior.
    for i in range (n+2):
        k = i
        a[k,k] = 1.0
        b[k] = bound (x[i],y[j],t)

    i = n+1                       # Frontera derecha.
    for j in range (m+2):
        k = j * (n+2) + i
        a[k,k] = 1.0
        b[k] = bound (x[i],y[j],t)
    
    j = m+1                       # Frontera superior.
    for i in range (n+2):
        k = j * (n+2) + i
        a[k,k] = 1.0
        b[k] = bound (x[i],y[j],t)
    
    i = 0                         # Frontera izquierda.
    for j in range (m+2):
        k = j * (n+2) + i
        a[k,k] = 1.0
        b[k] = bound (x[i],y[j],t)
    
    # -------------------------------------------------------------------
    # Resolución del sistema de ecuaciones.

    u = np.linalg.solve (a,b)

    plt.clf ()
#    titulo = 't = %5.2f' % t
#    graf01.graf (xx,yy,u,1,titulo)
#    print(len(xx))
#    print(len(yy))
#    print(len(u))

    print(u)

#    surf = ax.plot_surface(X=xx, Y=yy, Z=u, rstride=1, cstride=1, 
#cmap=cm.coolwarm, linewidth=0, antialiased=False, vmin=min(u), vmax=max(u))
#    ax.set_zlim(min(u), max(u))
    ax = fig.gca(projection='3d')

    Zu=np.reshape(u, (12,12))

    print(Zu)

#    surf = ax.plot_wireframe(X=xx, Y=yy, Z=Zu)
    ax.plot_surface(xx, yy, Zu, rstride=1, cstride=1, color='b')
    


#    ax.zaxis.set_major_locator(LinearLocator(10))
#    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=False)

#    fig.colorbar(surf)
#    plt.draw ()
    plt.draw()
    plt.show()
#    input ('?')
#    time.sleep (3.0)
    plt.pause(3)


    # Actualización.

    uold = u[:]

# -------------------------------------------------------------------
# Representación gráfica.

#plt.show ()

# -------------------------------------------------------------------

