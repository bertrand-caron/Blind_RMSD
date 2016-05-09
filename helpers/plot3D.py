#!/usr/bin/env python
import numpy.linalg
import pylab
import mpl_toolkits.mplot3d 

#draw a vector
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import math
import numpy as np


_fig = None
_ax = None

def plotPoint(point, color, marker, label):
    x, y, z = point
    ax = getFig()
    ax.scatter3D([x], [y], [z], color=color, marker=marker, s=100., label=label)

def plotPoints(points, color, marker, label):
    x, y, z = np.array(points).T
    ax = getFig()
    ax.scatter3D(x, y, z, color=color, marker=marker, s=100., label=label)

def plotVect(R, V):
    V = V/numpy.linalg.norm(V)
    V = R+V
    
    ax = getFig()
    
    a = Arrow3D([R[0], V[0]], [R[1], V[1]], [R[2], V[2]], mutation_scale=20, lw=2, arrowstyle="-|>", color="r")
    ax.add_artist(a)


def plotLine(line):
    ax = getFig()
    x, y, z = zip(*line)
    ax.plot(x,y,z)

def getFig():
    global _fig, _ax
    if _fig is not None and _ax is not None:
        return _ax
    else:
        _fig = pylab.figure(1)
        _fig.hold(True)
        _ax = mpl_toolkits.mplot3d.Axes3D(_fig)
        return _ax

def showGraph():
    global _fix, _ax
    _ax.legend(loc = 'upper left', numpoints=1, scatterpoints=1, frameon = False)
    # plot axes
    xlim = _ax.get_xlim3d()
    ylim = _ax.get_ylim3d()
    zlim = _ax.get_zlim3d()
    _ax.plot3D(xlim,[0,0],[0,0], color ='k', linestyle="--")
    _ax.plot3D([0,0],ylim,[0,0], color ='k', linestyle="--")
    _ax.plot3D([0,0],[0,0],zlim, color ='k', linestyle="--")
    _ax.set_xlabel("x")
    _ax.set_ylabel("y")
    _ax.set_zlabel("z")
    
    pylab.show()
    _fix, _ax = None, None

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

if __name__=="__main__":
    coords_ref = [ [1.,1.,0.], [1.,-1.,0.], [-1.,1.,0.], [-1.,-1.,0.]]
    coords_target = [[0.,math.sqrt(2),0.], [0.,-math.sqrt(2),0.], [math.sqrt(2),0.,0.], [-math.sqrt(2),0.,0.]]
    
    plotPoints(coords_ref, "r", "o", "ref")
    plotPoints(coords_target, "b", "^", "target")

    showGraph()
