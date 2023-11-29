from fenics import *
from mshr import *
from pylab import plt
import numpy as np

scale = 1/100

H  = 207.0 * scale
H1 = 250.0 * scale
H2 = 368.0 * scale

r  = 16.0 * scale
r0 = r/4
r1 = H/8

domain1 = Circle(Point(0,H),r,16)
domain2 = Rectangle(Point(-r/3, H),Point(r/3,H1))
domain3 = Rectangle(Point(-r/7, H),Point(r/7,H2))

ys = np.linspace(0,H,20)

poly = []

for y in ys:
    poly.append(Point(+(r0+(r-r0)*exp(-y/r1)),y))
for y in np.flip(ys,0):
    poly.append(Point(-(r0+(r-r0)*exp(-y/r1)),y))
    
base   = Polygon(poly)
domain = domain1 + domain2 + domain3 + base
mesh   = generate_mesh(domain,64)

mesh_file = File("tvtower.xml")
mesh_file << mesh
