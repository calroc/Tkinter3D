#!/usr/bin/env python
#
#    Copyright 2011 Simon Forman
#
#    This file is part of Tkinter3D.
#
#    Tkinter3D is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Tkinter3D is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Tkinter3D.  If not, see <http://www.gnu.org/licenses/>.
#
from Tkinter import Tk
from canvas3d import Canvas3D, Frame3D, roty
from scene import TkinterCanvasThing3D as dot


root = Tk()
c = Canvas3D(root)


# Zoom out a little.
c.frame.T.z += -300.0


# Create a dot at world origin (0, 0, 0).
origin = dot(c)
c.frame.things.append(origin)


# Make a frame for some objects to be in.
cube_frame = Frame3D()
c.frame.subframes.append(cube_frame)

# Add a Z component to the cube frame's translation to offset the cube
# from the world frame.
cube_frame.T.z += 300.0


# Make a cube.
F = -1, 1
corners = (
    (x, y, z)
    for x in F
        for y in F
            for z in F
    )
for x, y, z in corners:
    d = dot(c, 100.0 * x, 100.0 * y, 100.0 * z)
    cube_frame.things.append(d)


# Apply a rotation repeatedly to our cube.    
N = roty(360/30/3) # 4 degrees.
def delta():
    cube_frame.RM *= N
    c.after(60, delta)

# Start everything running.
delta()
c.startUpdating()
root.mainloop()
