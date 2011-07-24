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
Classes and functions for doing 3D math with quaternions.
'''
from math import sin, cos, pi, sqrt, acos, atan2, asin
from operator import add, sub, mul, div


tolerance = 0.00001 # A value "close enough" to zero.
scalar_types = int, float, long
radians = 180 / pi


def r2d(n):
    return n * radians


class Quaternion:

    def __init__(self, n=0.0, x=0.0, y=0.0, z=0.0, euler=False):
        self.n = n

        if isinstance(x, Vector):
            assert not euler #if euler, x,y,z should be angles
            self.v = x

        elif euler:
            #create Q from Euler angles
            x, y, z = (
                n / radians / 2.0
                for n in (x, y, z)
                )

            Q = (
                Quaternion(cos(x), sin(x))
                *
                Quaternion(cos(y), 0.0, sin(y))
                *
                Quaternion(cos(z), 0.0, 0.0, sin(z))
                )
            self.n = Q.n
            self.v = Q.v

        else:
            self.v = Vector(x, y, z)

    def __repr__(self):
        return "<Quaternion (%g, %s>" % (self.n, self.v.__repr__())

    def magnitude(self):
        v = self.v
        return sqrt(sum((
            self.n * self.n,
            v.x * v.x,
            v.y * v.y,
            v.z * v.z
            )))

    def getVector(self):
        return Vector(self.v.x, self.v.y, self.v.z)

    def getScalar(self):
        return self.n

    def __add__(self, p):
        if not isinstance(p, Quaternion):
            raise TypeError
        return Quaternion(self.n + p.n, self.v + p.v)

    def __sub__(self, p):
        if not isinstance(p, Quaternion):
            raise TypeError
        return Quaternion(self.n-p.n, self.v-p.v)

    def __mul__(self, p):
        if isinstance(p, scalar_types):
            return Quaternion(self.n * p, self.v * p)

        elif isinstance(p, Quaternion):
            return Quaternion(
                self.n * p.n - self.v * p.v,  # the real part "n"
                (
                    p.v * self.n
                    + self.v * p.n
                    + self.v ^ p.v
                    )
                )

        elif isinstance(p, Vector):
            return self * Quaternion(0.0, p)

        else:
            raise TypeError, p

    def __div__(self, p):
        if not isinstance(p, scalar_types):
            raise TypeError
        return Quaternion(self.n / p, self.v / p)

    def __invert__(self):
        """Return the conjugate itself."""
        return Quaternion(self.n, ~self.v)

    def getAngle(self):
        """Angle of rotation about axis."""
        return 2 * acos(self.n)

    def getAxis(self):
        m = self.v.magnitude()
        if m <= tolerance:
            return Vector()
        else:
            return self.v / m

    def Rotate(q, p):
        if not isinstance(p, Quaternion):
            raise TypeError
        return q * p * (~q)

    def vectorRotate(q, v):
        if not isinstance(v, Vector):
            raise TypeError
        t = q * v * (~q)
        return t.getVector()

    def getEulerAngles(q):
        q00 = q.n * q.n
        q11 = q.v.x * q.v.x
        q22 = q.v.y * q.v.y
        q33 = q.v.z * q.v.z

        r11 = q00 + q11 - q22 - q33
        r21 = 2 * (q.v.x * q.v.y + q.n * q.v.z)
        r31 = 2 * (q.v.x * q.v.z + q.n * q.v.y)
        r32 = 2 * (q.v.y * q.v.z + q.n * q.v.x)
        r33 = q00 - q11 - q22 + q33

        tmp = abs(r31)
        if tmp > 0.999999:
            r12 = 2 * (q.v.x * q.v.y - q.n * q.v.z)
            r13 = 2 * (q.v.x * q.v.z + q.n * q.v.y)

            return Vector(
                0.0,
                r2d(-(pi / 2) * (r31 / tmp)),
                r2d(atan2(-r12, -r31*r13))
                )
        return Vector(
            r2d(atan2(r32, r33)),
            r2d(asin(-r31)),
            r2d(atan2(r21, r11))
            )


class Vector:
    """
    This defines a vector of three coordinates.
    """

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "<Vector (%g,%g,%g)>" % (self.x, self.y, self.z)

    def magnitude(self):
        return sqrt(
            self.x * self.x +
            self.y * self.y +
            self.z * self.z
            )

    def normalize(self):
        m = self.magnitude()
        if m < tolerance:
            m = 1
        xyz = [self.x / m, self.y / m, self.z / m]
        for i in (0, 1, 2):
            if abs(xyz[i]) < tolerance:
                xyz[i] = 0.0
        return Vector(*xyz)

    def reverse(self):
        return Vector(-self.x, -self.y, -self.z)

    def __add__(self, r_v):
        if not isinstance(r_v, Vector):
            raise TypeError
        return Vector(self.x + r_v.x, self.y + r_v.y, self.z + r_v.z)

    def __sub__(self, r_v):
        if not isinstance(r_v, Vector):
            raise TypeError
        return Vector(self.x - r_v.x, self.y - r_v.y, self.z - r_v.z)

    def __mul__(self, r_n):
        if isinstance(r_n, Vector):
            #Dot Product
            return sum((
                self.x * r_n.x,
                self.y * r_n.y,
                self.z * r_n.z
                ))

        elif isinstance(r_n, Matrix):
            #Matrix multiplication
            return r_n.__mul__(self)

        elif isinstance(r_n, scalar_types):
            #Scalar multiplication
            return Vector(self.x * r_n, self.y * r_n, self.z * r_n)

        else:
            raise TypeError, r_n

    def __div__(self, r_n):
        if not isinstance(r_n, scalar_types):
            raise TypeError
        return Vector(self.x / r_n, self.y / r_n, self.z / r_n)

    def __neg__(self):
        return self.reverse()

    def __invert__(self):
        return self.reverse()

    def __xor__(self, r_v):
        """Note: NOT eXclusive-OR, Vector Cross Product."""
        if not isinstance(r_v, Vector):
            raise TypeError
        return Vector(
            self.y * r_v.z - self.z * r_v.y,
           -self.x * r_v.z + self.z * r_v.x,
            self.x * r_v.y - self.y * r_v.x
            )

    def triplescalar(self, v, w):
        if not (isinstance(v, Vector)
                and isinstance(w, Vector)):
            raise TypeError
        res = self * (v ^ w)
        if abs(res) < tolerance:
            res = 0.0
        return  res


class Matrix:

    make_local_variables = """\
e11=self.e11; e12=self.e12; e13=self.e13
e21=self.e21; e22=self.e22; e23=self.e23
e31=self.e31; e32=self.e32; e33=self.e33
"""

    def __init__(self,
                 e11=0.0,  e12=0.0,  e13=0.0,
                 e21=0.0,  e22=0.0,  e23=0.0,
                 e31=0.0,  e32=0.0,  e33=0.0
                 ):
        self.e11 = e11
        self.e12 = e12
        self.e13 = e13
        self.e21 = e21
        self.e22 = e22
        self.e23 = e23
        self.e31 = e31
        self.e32 = e32
        self.e33 = e33

    def __repr__(self):
        exec self.make_local_variables
        return "<Matrix (\n%f, %f, %f,\n%f, %f, %f,\n%f, %f, %f\n)>" % (
            e11,  e12,  e13,
            e21,  e22,  e23,
            e31,  e32,  e33
            )

    def determinant(self):
        exec self.make_local_variables
        return e11*e22*e33 - \
               e11*e32*e23 + \
               e21*e32*e13 - \
               e21*e12*e33 + \
               e31*e12*e23 - \
               e31*e22*e13

    def transpose(self):
        exec self.make_local_variables
        return Matrix(
            e11, e21, e31,
            e12, e22, e32,
            e13, e23, e33
            )

    def __invert__(self):
        """~M"""
        exec self.make_local_variables
        d = self.determinant()
        if d == 0.0:
            d = 1.0
        return Matrix(*[
            n / d
            for n in (
                e22*e33 - e23*e32, -(e12*e33 - e13*e32),  e12*e23 - e13*e22,
              -(e21*e33 - e23*e31),  e11*e33 - e13*e31, -(e11*e23 - e13*e21),
                e21*e32 - e22*e31, -(e11*e32 - e12*e31),  e11*e22 - e12*e21
                )
            ])
        # N.B. I had to add the '*[]' notation in the above call to
        # Matrix() because of a weird interaction with the exec:
        #
        # $ python math3d.py 
        #     File "math3d.py", line 298
        #       exec self.make_local_variables
        #   SyntaxError: unqualified exec is not allowed in function
        #   '__invert__' it contains a nested function with free variables
        #

    def __op(n, m, op): #Note the arg traditionally called 'self' is here called 'n'
        """This function applies the op 'op' to each element pair of two matricies."""
        ne = [
            n.e11,  n.e12,  n.e13,
            n.e21,  n.e22,  n.e23,
            n.e31,  n.e32,  n.e33
            ]
        me = [
            m.e11,  m.e12,  m.e13,
            m.e21,  m.e22,  m.e23,
            m.e31,  m.e32,  m.e33
            ]
        return Matrix(*map(op, ne, me))

    def __add__(self, r_v):
        if isinstance(r_v, scalar_types):
            # For scalar addition we convert the scalar value into a
            # matrix containing the scalar value in every element and
            # then add it, below.
            r_v = Matrix(*((r_v,) * 9))

        elif not isinstance(r_v, Matrix):
            raise TypeError, r_v

        return self.__op(r_v, add)

    def __sub__(self, r_v):
        if isinstance(r_v, scalar_types):
            r_v = Matrix(*((r_v,) * 9))

        elif not isinstance(r_v, Matrix):
            raise TypeError, r_v

        return self.__op(r_v, sub)

    def __mul__(self, r_n):
        if isinstance(r_n, Matrix):
            exec self.make_local_variables
            me11=r_n.e11; me12=r_n.e12; me13=r_n.e13
            me21=r_n.e21; me22=r_n.e22; me23=r_n.e23
            me31=r_n.e31; me32=r_n.e32; me33=r_n.e33

            return Matrix(
                e11*me11 + e12*me21 + e13*me31,
                e11*me12 + e12*me22 + e13*me32,
                e11*me13 + e12*me23 + e13*me33,

                e21*me11 + e22*me21 + e23*me31,
                e21*me12 + e22*me22 + e23*me32,
                e21*me13 + e22*me23 + e23*me33,

                e31*me11 + e32*me21 + e33*me31,
                e31*me12 + e32*me22 + e33*me32,
                e31*me13 + e32*me23 + e33*me33,
                )

        elif isinstance(r_n, Vector):
            exec self.make_local_variables
            return Vector(
                e11*r_n.x + e12*r_n.y + e13*r_n.z,
                e21*r_n.x + e22*r_n.y + e23*r_n.z,
                e31*r_n.x + e32*r_n.y + e33*r_n.z
                )

        elif isinstance(r_n, scalar_types):
            tmp = Matrix(*((r_n,) * 9))
            return self.__op(tmp, mul)

        else:
            raise TypeError, r_n

    def __div__(self, r_n):
        if not isinstance(r_n, scalar_types):
            raise TypeError
        return self.__op(Matrix(*((r_n,) * 9)), div)


def planeNormalAndDistance(v0, v1, v2):
    pN = (v2 - v0) ^ (v1 - v0)
    pN = pN.normalize()
    d  = v0 * pN
    return pN, d


def Frustum2Planes(f):
    v0, v1, v2, v3, v4, v5, v6, v7 = f
    lpN, ld = planeNormalAndDistance(v0, v4, v3)
    rpN, rd = planeNormalAndDistance(v6, v2, v5)
    tpN, td = planeNormalAndDistance(v3, v7, v1)
    bpN, bd = planeNormalAndDistance(v2, v6, v0)
    return (lpN, ld), (rpN, rd), (tpN, td), (bpN, bd)


def distanceVector2Plane(pN, pd, v):
    return pN * v - pd


def distanceVector2PlanePlus(pN, pd, v):
    return pN.x*v.x + pN.y*v.y + pN.z*v.z - pd



def rotx(a):
    a = a / radians
    return Matrix(
        1,      0,       0,
        0, cos(a), -sin(a),
        0, sin(a),  cos(a)
        )


def roty(a):
    a = a / radians
    return Matrix(
        cos(a), 0, sin(a),
        0,      1,      0,
       -sin(a), 0, cos(a)
        )


def rotz(a):
    """Clockwise as facing z+ """
    a = a / radians
    return Matrix(
        cos(a), -sin(a), 0,
        sin(a),  cos(a), 0,
        0,            0, 1
        )
