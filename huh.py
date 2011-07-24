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
from math3d import (
    scalar_types,
    Vector,
    Matrix,
    rotx, roty, rotz,
    )


class Frame3D:
    '''
    A frame, has scale, translation, and rotation.
    '''
    def __init__(self,
        x=0.0,
        y=0.0,
        z=0.0,
        scale=1.0,
        translation=None,
        rotation=None,
        ):
        for n in (x, y, z, scale):
            if not isinstance(n, scalar_types):
                raise TypeError, n

        self.s = scale

        # Translation.
        if isinstance(translation, Vector):
            self.T = translation

        elif isinstance(translation, (tuple, list)):
            if not len(translation) == 3:
                raise ValueError, translation
            for t in translation:
                if not isinstance(t, scalar_types):
                    raise TypeError, t
            self.T = Vector(*translation)

        else:
            self.T = Vector()

        # Rotation.
        if isinstance(rotation, Matrix):
            self.RM = rotation

        elif isinstance(rotation, (tuple, list)):
            if not len(rotation) == 3:
                raise ValueError, rotation
            for a in rotation:
                if not isinstance(a, scalar_types):
                    raise TypeError, a
            self.RM = (
                rotx(rotation[0])
                * roty(rotation[1])
                * rotz(rotation[2])
                )

        else:
            self.RM = rotx(0) * roty(0) * rotz(0)

        self.vector_set = [Vector(x, y, z)]
        self.items = {} # Canvas item id for our vector/dot.
        self._W = 3 # width in pixels of canvas dot.
        self.canvas.render_list.append(self)

    def H(self, v):
        '''
        Apply scale, translation, and rotation defined here to a Vector.
        '''
##        return v + self.T
    # adding two vectors is an order of magnitude faster than
    # multiplying a vector by a matrix.
        return (v * self.s - self.T) * self.RM


    def render(self):
        '''
        Manage this dot on its Canvas.
        '''

        # Get our width and vector.
        W = self._W
        v = self.vector_set[0]

        # Compute our transform.
        wv = self.H(v) 

        if self.canvas.visible(wv):

            # Get our x, y "screen" coords and our canvas item id.
            x, y = self.canvas.vectorToScreen(wv)
            item = self.items.get(v, None)

            if item:
                # Move our item if we have one.
                self.canvas.coords(item, x - W, y - W, x + W, y + W)

            else:
                # Make an oval item and remember its id.
                self.items[v] = self.canvas.create_oval(
                    x - W, y - W, x + W, y + W,
                    fill='green'
                    )
        else:
            # We're not visible, remove our oval if any.
            item = self.items.get(v, None)
            if item:
                self.canvas.delete(item)
                self.items[v] = None


