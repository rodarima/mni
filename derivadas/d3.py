# encoding: utf-8

import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.Qt import QtCore, QtGui

tfin = 1
nt = 200
n = 20
m = 20

xini = 0.0
xfin = 1.0
yini = 0.0
yfin = 1.0

nodos = (n + 2) * (m + 2)

x = np.linspace (xini, xfin, n+2)
y = np.linspace (yini, yfin, m+2)
incx = x[1] - x[0]
incy = y[1] - y[0]
inct = tfin / nt

xx,yy = np.meshgrid(x,y)

def frontera(x, y, t):
#    if(x==0 and t>tfin/2): return 1.0
#    if(x==0 and t<tfin/2): return 1.0
    if(x==0):
#        if(t<tfin/4): return t/(tfin/4)
#        if(t<tfin/2): return 1.0
#        if(t<tfin*3/4): return (tfin*3/4-t)/(tfin/4)
        if(t<tfin/2): return 1.0
#    if(x==0.0): return 1.0
    return 0.0

def flujo(x, y, t):
    if(t > tfin/4.0): return -1.0
    return 0.0

def inicial (x,y):
#    return x * y * (1.0 - x) * (1.0 - y) * 16
#    if(x<0.1): return 1.0
#    return x
#    return 1.0
    return 0.0

def f(x,y,t):
#    return x * y * (1.0 - x) * (1.0 - y) * 16
#    return x
    return 1.0

# Matriz de coeficientes de la solución inicial. U0
S = np.zeros(nodos, 'f')

#a = 0.993
a = 1.0
cx2 = a/incx**2.
cy2 = a/incy**2.
cuk = 1./inct + 2.*a*cx2 + 2.*a*cy2

U = np.zeros([nodos, nodos], 'f')
B = np.zeros(nodos, 'f')

S = np.zeros(nodos, 'f')
Z = np.zeros([nt, n+2, m+2], 'f')

for i in range (n+2):
    for j in range (m+2):
        k = j * (n+2) + i
        S[k] = inicial(x[i], y[j])

for it in np.arange(nt):
    t = (it + 1) * inct
    per = it * 100 / nt
    print(str(per) + "%")
    for i in np.arange(1,n+1):
        for j in np.arange(1,m+1):

            k = j*(n+2) + i

            B[k] = 1.0/inct * S[k] + f(x[i], x[j], t)
            U[k,k] = cuk
            U[k,k+1] = -cx2
            U[k,k-1] = -cx2
            U[k,k+(m+2)] = -cy2
            U[k,k-(m+2)] = -cy2

    for i in np.arange(n+2):
        j = 0
        k = i
        U[k,k] = 1.0
        B[k] = frontera(x[i], y[j], t)

        j = m+1
        k = j * (n+2) + i
        U[k,k] = 1.0
        B[k] = frontera(x[i], y[j], t)
    
    for j in np.arange(m+2):
        i = 0
        k = j * (n+2) + i
        U[k,k] = 1.0
        B[k] = frontera(x[i], y[j], t)
        #Frontera con flujo
        i = n+1
        k = j * (n+2) + i
        U[k,k] = 1.0 / incx
        U[k,k-1] = -1.0 / incx
        B[k] = flujo(x[i], y[j], t)
    
    S = np.linalg.solve(U,B)
    S = np.array(S)


    Z[it] = np.reshape(S, (n+2,m+2))

app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: GLVolumeItem')
#gls = gl.GLSurfacePlotItem(x=x, y=y, z=Z[2], shader='shaded')
g = gl.GLGridItem()
g.scale(.1,.1,.1)
g.setDepthValue(10)
w.addItem(g)
#gls = gl.GLSurfacePlotItem(x=x, y=y, z=Z[0], shader='heightColor')
#gls.shader()['colorMap'] = np.array([0.2, 2, 0.5, 0.2, 1, 1, 0.2, 0, 2])
gls = gl.GLSurfacePlotItem(x=x, y=y, z=Z[0], shader='shaded', 
        glOptions='opaque')
gls.scale(1,1,1)
gls.translate(-.5, -.5, 0.)
w.addItem(gls)

index = 0
def update():
    global index, Z, nt
    index+=1
    gls.setData(z=Z[index%nt])
    
timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(50)

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()


