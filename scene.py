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
        scale=1.0,
        translation=None,
        rotation=None,
        ):

        if not isinstance(scale, scalar_types):
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

        self.things = []
        self.subframes = []

    def transform(self, vector):
        '''
        Apply scale, translation, and rotation defined here to a Vector.
        '''
        if self.s != 1:
            vector = vector * self.s
        return (vector - self.T) * self.RM

    def yieldTransformed(self):
        for subframe in self.subframes:
            for thing, vector in subframe.yieldTransformed():
                yield thing, self.transform(vector)
        for thing in self.things:
            yield thing, self.transform(thing.getVector())


class Thing3D:
    '''
    A thing that has a location in three-space.
    '''
    def __init__(self, x=0.0, y=0.0, z=0.0):
        for n in (x, y, z):
            if not isinstance(n, scalar_types):
                raise TypeError, n

        self.vector = Vector(x, y, z)

    def getVector(self):
        return self.vector

    def hide(self):
        pass

    def render(self, x, y, z):
        pass


class TkinterCanvasThing3D(Thing3D):

    def __init__(self, canvas, x=0.0, y=0.0, z=0.0, **settings):
        Thing3D.__init__(self, x, y, z)

        self.canvas = canvas

        # Width in pixels of canvas dot.
        self._width = settings.pop('width', None) or 3

        self.settings = settings
        self.settings.setdefault('fill', 'green')

        # Canvas item id for our vector/dot.
        self.item = None

    def hide(self):
        if self.item:
            self.canvas.delete(self.item)
            self.item = None

    def render(self, x, y, z):
        '''
        Manage this dot on its Canvas.
        '''
        # Get our width and vector.
        W = self._width

        # Move our item if we have one.
        if self.item:
            self.canvas.coords(self.item, x - W, y - W, x + W, y + W)

        # Make an oval item and remember its id.
        else:
            self.item = self.canvas.create_oval(
                x - W, y - W, x + W, y + W,
                **self.settings
                )


