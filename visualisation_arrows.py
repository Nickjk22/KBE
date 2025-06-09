#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020 ParaPy Holding B.V.
#
# This file is subject to the terms and conditions defined in
# the license agreement that you have received with this source code
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR
# PURPOSE.

from parapy.core import Attribute, Input, Part, derived
from parapy.geom import (
    Compound, Cone, Cylinder, GeomBase, Orientation, Position, translate)
from parapy.geom.generic.positioning import orthogonal_vector


class ArrowVisualisationMixin(GeomBase):
    """Mixin to add a Arrow visualization to a Load or Support."""

    #: The dimension for visualization. Corresponds to length of an arrow.
    #: :type: float
    v_dim = Input(1.0)

    #: Color of arrows.
    #: :type: str | (float, float, float)
    color = Input("red")

    @Attribute
    def arrow_anchor_points(self):
        """:rtype: list[parapy.geom.Point]"""
        raise NotImplementedError()

    @Attribute
    def _visual_objects(self):
        """:rtype: list[parapy.geom.occ.drawable.DrawableShape]"""
        raise NotImplementedError()

    @Part
    def visual(self):
        return Compound(self._visual_objects)


class VisualizationBase(Compound):
    #: :type: parapy.geom.Point
    point = Input()

    #: :type: parapy.geom.Vector
    direction = Input()

    #: pointing away from :attr:`point` or pointing towards it?
    #: :type: bool
    outward = Input(True)


from parapy.core import *
from parapy.geom import *
from math import isclose

class LiftArrowArray(GeomBase):
    __initargs__ = ["points_list", "lift_forces"]

    points_list = Input()
    lift_forces = Input()

    nb_of_heads = Input(1)

    @Attribute
    def max_force(self):
        """Maximale absolute kracht om te normaliseren."""
        if not self.lift_forces:
            return 1.0  # voorkomt deling door nul
        return max(abs(f) for f in self.lift_forces if not isclose(f, 0.0))


    @Part
    def arrows(self):
        return CylinderArrow(
            quantify=len(self.points_list),
            point=self.points_list[child.index],
            direction=Vector(0, 0, 1) if self.lift_forces[child.index] >= 0 else Vector(0, 0, -1),
            head_length=abs(self.lift_forces[child.index]) / self.max_force,
            nb_of_heads=self.nb_of_heads,
        )

class CylinderArrow(VisualizationBase):
    __initargs__ = ["point", "direction", "head_length"]

    point = Input()
    direction = Input()
    head_length = Input(1.0)
    nb_of_heads = Input(1)
    outward = Input(True)

    @Attribute
    def base_length(self):
        return self.head_length * 3

    @Attribute
    def head_radius(self):
        return 0.5 * self.head_length

    @Attribute
    def base_radius(self):
        return 0.2 * self.head_radius

    @Attribute
    def built_from(self):

        pt = self.point
        v = self.direction.normalize
        n = self.nb_of_heads
        lh = self.head_length
        lbase = self.base_length
        rbase = self.base_radius

        if not self.outward:
            pt = pt - v * (lbase + n * lh)

        w = orthogonal_vector(v)
        pos = Position(pt, Orientation(x=w, z=v))

        lst = [Cylinder(radius=rbase, height=lbase, position=pos)]
        rhead = self.head_radius

        for i in range(n):
            pos_ = translate(pos, 'z', lbase + i * lh)
            cone = Cone(radius1=rhead, radius2=0, height=lh, position=pos_)
            lst.append(cone)
        return lst
