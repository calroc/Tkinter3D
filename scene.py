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
    rotx,
    roty,
    rotz,
    )


class Frame3D:
    '''
    A frame of reference for applying a scale, translation, and rotation
    to a collection of :class:`Thing3D` objects. It can also have
    descendant :class:`Frame3D` objects that contain their own
    :class:`Thing3D` objects.

    :param scale: A scaling factor.
    :type scale: One of :data:`math3d.scalar_types`

    :param translation: Specifies XYZ translation.
    :type translation: :class:`math3d.Vector` or any iterable of any
       three :data:`math3d.scalar_types`

    :param rotation: Specifies a rotation around each of the X, Y, and Z
       axises. If given as discrete values (i.e. not a
       :class:`math3d.Matrix`) they are interpreted as degrees.
    :type rotation: :class:`math3d.Matrix` or any iterable of any
       three :data:`math3d.scalar_types`

    :raises: :exc:`TypeError` if the constraints on parameters are not
       met.

    Each of the above three parameters becomes an instance attribute of
    the resulting :class:`Frame3D` object that you can manipulate at
    runtime to affect subsequent calls to :meth:`yieldTransformed`.


    .. attribute:: s

        (Instance attribute.) Scale factor applied by
        :meth:`yieldTransformed`.

    .. attribute:: T

        (Instance attribute.) :class:`math3d.Vector` Translation applied
        by :meth:`yieldTransformed`.

    .. attribute:: RM

        (Instance attribute.) :class:`math3d.Matrix` Rotation applied by
        :meth:`yieldTransformed`.

    .. attribute:: things

        (Instance attribute.) A list of :class:`Thing3D` objects.  Append
        :class:`Thing3D` objects to this list to have them yielded by
        :meth:`yieldTransformed`.

    .. attribute:: subframes

        (Instance attribute.) A list of :class:`Frame3D` objects.  Append
        :class:`Frame3D` objects to this list to have their
        :class:`Thing3D` objects yielded by :meth:`yieldTransformed`.

    '''

    def __init__(
        self,
        scale=1.0,
        translation=None,
        rotation=None,
        ):

        self.things = []
        self.subframes = []

        if not isinstance(scale, scalar_types):
            raise TypeError, scale
        self.s = scale

        # Translation.
        if isinstance(translation, Vector):
            self.T = translation

        elif translation is not None:

            translation = tuple(translation)

            if not (
                len(translation) == 3
                and all(isinstance(t, scalar_types) for t in translation)
                ):
                raise TypeError, translation

            self.T = Vector(*translation)

        else:
            self.T = Vector()

        # Rotation.
        if isinstance(rotation, Matrix):
            self.RM = rotation

        elif rotation is not None:

            rotation = tuple(rotation)

            if not (
                len(rotation) == 3
                and all(isinstance(a, scalar_types) for a in rotation)
                ):
                raise TypeError, rotation

            ax, ay, az = rotation
            self.RM = rotx(ax) * roty(ay) * rotz(az)

        else:
            self.RM = rotx(0) * roty(0) * rotz(0)

    def transform(self, vector):
        '''
        Return a Vector that is the result of appling the scale,
        translation, and rotation defined by this :class:`Frame3D` to the
        given input Vector.

        :param vector: A Vector representing a 3D point.
        :type vector: :class:`math3d.Vector`
        :rtype: :class:`math3d.Vector`

        '''
        if self.s != 1:
            vector = vector * self.s
        return (vector - self.T) * self.RM

    def yieldTransformed(self):
        '''
        Yield each :class:`Thing3D` in this :class:`Frame3D` and its
        descendant :class:`Frame3D` objects, modified by its scale,
        translation, and rotation.
        This is called recursively on descendant :class:`Frame3D` objects.

        :rtype: generator of two-tuples:
            (:class:`Thing3D`, transformed :class:`math3d.Vector`)

        '''
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
        self._width = (settings.pop('width', None) or 3) / 2

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

        W = self._width

        top, left, bottom, right = x - W, y - W, x + W, y + W

        # Move our item if we have one.
        if self.item:
            self.canvas.coords(self.item, top, left, bottom, right)

        # Make an oval item and remember its id.
        else:
            self.item = self.canvas.create_oval(
                top,
                left,
                bottom,
                right,
                **self.settings
                )
