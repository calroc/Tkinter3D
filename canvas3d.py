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
'''
A Docstring for canvas3d module.
'''
from Tkinter import Tk, Canvas
from math3d import (
    Vector,
    planeNormalAndDistance,
    rotx,
    roty,
    rotz,
    radians,
    )
from math import tan, pi
from scene import Frame3D


def G(L, a):
    return L / tan(a / radians) #the arc of the viewport


class Canvas3D(Canvas):
    '''
    A Tkinter Canvas adapted to serve as a 3D viewer.
    '''

    def __init__(self, master=None, fps=20, *args, **kw):
        Canvas.__init__(self, *((master,) + args), **kw)

        self.frame = Frame3D()
        self.items = {}

        self.w = float(self['width'])
        self.h = float(self['height'])

        self.axes_on = True
        self.xyaxes()

        self.bind('<Configure>', self.resize_callback)

        self.fps_ms = 1000 / fps # Frames per second in milliseconds.

        self.focal_distance = 2.2
        self.x_arc_of_view = 73.0
        self.y_arc_of_view = 73.0

        self.f = G(self.w, self.x_arc_of_view)
        self.u = G(self.h, self.y_arc_of_view)

        self._frustum = Frustum(
            self.w, self.h,
            self.x_arc_of_view, self.y_arc_of_view,
            self.focal_distance,
            )

    def visible(self, v):
        '''
        Return a Boolean indicating if the given Vector is within the
        view.

        :param v: A Vector representing a 3D point.
        :type v: :class:`math3d.Vector`
        :rtype: Boolean

        '''
        return self._frustum.visible(v)

    def vectorToScreen(self, v):
        z = v.z
        if z == 0.0:
            z = 1 # z = tolerance?
        FD = self.focal_distance
        w = self.w / 2
        f = self.f 
        h = self.h / 2
        u = self.u

        return (
            int(w + v.x * f / (z + f * FD)),
            int(h - v.y * u / (z + u * FD))
            )

    def xyaxes(self):
        if not self.axes_on:
            return

        w, h = int(self['width']), int(self['height'])

        yaxis = self.items.get('yaxis', None)
        if yaxis:
            self.coords(yaxis, w / 2, 0, w / 2, h)
        else:
            self.items['yaxis'] = self.create_line(w / 2, 0, w / 2, h)

        xaxis = self.items.get('xaxis', None)
        if xaxis:
            self.coords(xaxis, w, h / 2, 0, h / 2)
        else:
            self.items['xaxis'] = self.create_line(w, h / 2, 0, h / 2)

    def resize_callback(self, event):
        self.xyaxes()

    def startUpdating(self):
        self._updater()

    def _updater(self):
        for thing, vector in self.frame.yieldTransformed():

            if self._frustum.visible(vector):
                x, y = self.vectorToScreen(vector)
                thing.render(x, y, 1.0)

            else:
                thing.hide()

        self.update_idletasks()
        self.update()
        self.after(self.fps_ms, self._updater)


class Frustum:

    def __init__(
        self,
        width=320,
        height=240,
        x_arc=73.0,
        y_arc=73.0,
        focal_distance=2.2,
        depth=5000.0,
        ):

        self.w, self.h = width, height
        self.xarc, self.yarc = x_arc, y_arc
        self.FD = focal_distance
        self.depth = depth

        self.T = Vector()
        self.RM = rotx(0) * roty(0) * rotz(0)

        # Create a list of eight "blank" new Vectors. These will be the
        # corners of the frustum, used to determine the planes beyond
        # which a given point is outside the viewable space.
        self._frustum = [Vector() for _ in range(8)]

        self._reset()
        self._getFrustumPlanes()

    def _reset(self):
        self.f = G(self.w, self.xarc)
        self.u = G(self.h, self.yarc)
        self._frustum_dirty = True
        self._frustum_planes_dirty = True

    def transform(self, v):
        return (v - self.T) * self.RM

    def visible(self, v):
        """
        Is vector v in World coordinates within the view frustum?
        """
        (lpN, ld), (rpN, rd), (tpN, td), (bpN, bd) = self._getFrustumPlanes()
        vx, vy, vz = v.x, v.y, v.z
        accs = acl, acr, act, acb = -ld, -rd, -td, -bd

        if vx != 0.0:
            acl += vx * lpN.x
            acr += vx * rpN.x
            act += vx * tpN.x
            acb += vx * bpN.x

        if vy != 0.0:
            acl += vy * lpN.y
            acr += vy * rpN.y
            act += vy * tpN.y
            acb += vy * bpN.y

        if vz != 0.0:
            if (
                (acl + vz * lpN.z <= 0.0)
                or (acr + vz * rpN.z <= 0.0)
                or (act + vz * tpN.z <= 0.0)
                ):
                return False
            return acb + vz * bpN.z > 0.0

        else:
            if acl <= 0.0: return False
            if acr <= 0.0: return False
            if act <= 0.0: return False
            if acb >  0.0: return True
            return False

    def _getFrustum(self):
        """
        This returns the coordinates of the view frustum in self/model
        coordinates.
        """
        if self._frustum_dirty:
            FD, w, h, f, u = self.FD, self.w, self.h, self.f, self.u
            v0, v1, v2, v3, v4, v5, v6, v7 = self._frustum
            z = min(f, u)

            xw =  w / 2.0 * (z + f * FD) / f
            x0 = -xw
            yh =  h / 2.0 * (z + u * FD) / u
            y0 = -yh

            v0.x, v0.y, v0.z = x0, y0, z
            v1.x, v1.y, v1.z = xw, yh, z
            v2.x, v2.y, v2.z = xw, y0, z
            v3.x, v3.y, v3.z = x0, yh, z

            z = self.depth
            xw =  w / 2.0 * (z + f * FD) / f
            x0 = -xw
            yh =  h / 2.0 * (z + u * FD) / u
            y0 = -yh

            v4.x, v4.y, v4.z = x0, y0, z
            v5.x, v5.y, v5.z = xw, yh, z
            v6.x, v6.y, v6.z = xw, y0, z
            v7.x, v7.y, v7.z = x0, yh, z

            self._frustum_dirty = False

        return self._frustum

    def _getFrustumPlanes(self):
        if self._frustum_planes_dirty:

            V = map(self.transform, self._getFrustum())

            self._frustum_planes = tuple(
                planeNormalAndDistance(V[i0], V[i1], V[i2])
                for i0, i1, i2 in (
                    (0, 4, 3),
                    (6, 2, 5),
                    (3, 7, 1),
                    (2, 6, 0),
                    )
                )
            self._frustum_planes_dirty = False

        return self._frustum_planes


if __name__ == '__main__':
    from scene import Thing3D as dot

    root = Tk()
    root.title("Rotating Cube Demo")
    c = Canvas3D(root)
    c.pack(expand=1, fill='both')

## Changing the frustum's "frame".
#    c._frustum.T.z += 100.0
#    c._frustum._reset()

    # Zoom out a little.
    c.frame.T.z += -200.0

    # Create a dot at world origin (0, 0, 0).
    origin = dot(c)
    c.frame.things.append(origin)

    # Make a cube.
    cube_frame = Frame3D()
    c.frame.subframes.append(cube_frame)
    coords = -1, 1
    for x in coords:
        for y in coords:
            for z in coords:
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
