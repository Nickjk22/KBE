#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2021 ParaPy Holding B.V.
#
# This file is subject to the terms and conditions defined in
# the license agreement that you have received with this source code
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR
# PURPOSE.

from typing import Optional, Sequence

from OCC.utils.top import create_toptools_listofshape
from OCC.wrapper.BOPAlgo import (
    BOPAlgo_CellsBuilder, BOPAlgo_GlueFull, BOPAlgo_GlueOff, BOPAlgo_GlueShift)

from parapy.core import Attribute, Input
from parapy.core.widgets import CheckBox, Dropdown
from parapy.geom.occ.brep import BRep, BRepBuilder
from parapy.geom.occ.compound import Compound_

__all__ = ["GeneralFuse"]

PY_TO_BOPAlgo_GlueEnum = {"off": BOPAlgo_GlueOff,
                          "shift": BOPAlgo_GlueShift,
                          "full": BOPAlgo_GlueFull}


class GeneralFuse(Compound_, BRepBuilder):
    """General Fuse Algorithm. Its main purpose is to build the split parts
    of the :attr:`shape_in` and :attr:`tools` from which the result of the
    operations is combined. Usage:

    >>> from parapy.geom import RectangularFace
    >>> s1 = RectangularFace(1, 1)
    >>> s2 = RectangularFace(1, 1, position=s1.position.translate(x=0.5, y=0.5))
    >>> obj = GeneralFuse(tools=[s1, s2])
    >>> len(obj.faces), len(obj.modified(s1)), len(obj.modified(s2))
    (3, 2, 2)

    The above 2 faces where creating 3 new faces, 2 modified from s1 and
    2 modified from s2, 1 face is common to both.

    By default all splits parts will become part of the result, however
    :attr:`to_keep` and :attr:`to_avoid` can influence this. If specified,
    to be taken into the result, a part must be IN for all shapes from the
    list :attr:`to_keep` and must be OUT of all shapes from the list
    :attr:`to_avoid`.

    >>> # Subtract(s1, s2)
    >>> obj.to_keep = [s1]
    >>> obj.to_avoid = None  # the default value, but shown for consistency
    >>> len(obj.faces), len(obj.modified(s1)), len(obj.modified(s2))
    (2, 2, 1)
    >>> # Subtract(s2, s1)
    >>> obj.to_keep = [s2]
    >>> obj.to_avoid = None  # the default value, but shown for consistency
    >>> len(obj.faces), len(obj.modified(s1)), len(obj.modified(s2))
    (2, 1, 2)
    >>> # Common(s1, s2)
    >>> obj.to_keep = [s1, s2]
    >>> obj.to_avoid = None  # the default value, but shown for consistency
    >>> len(obj.faces), len(obj.modified(s1)), len(obj.modified(s2))
    (1, 1, 1)
    >>> # keep s1, but avoid s2
    >>> obj.to_keep = [s1]
    >>> obj.to_avoid = [s2]
    >>> len(obj.faces), len(obj.modified(s1)), len(obj.modified(s2))
    (1, 1, 0)
    >>> # keep s2, but avoid s1
    >>> obj.to_keep = [s2]
    >>> obj.to_avoid = [s1]
    >>> len(obj.faces), len(obj.modified(s1)), len(obj.modified(s2))
    (1, 0, 1)

    The API envisions 2 usage scenarios:

    1. It you have some basis shape that you want to intersect with *tools* ::

       >>> obj = GeneralFuse(shape_in=S, tools=[T1, ..., TN])  # doctest: +SKIP

    2. All shapes are equal, just use ``tools``::

       >>> obj = GeneralFuse(tools=[S1, ..., SN])  # doctest: +SKIP
    """
    __initargs__ = ("shape_in", "tools")

    shape_in: Optional[BRep] = Input(None)
    tools: Sequence[BRep] = Input()
    glue = Input("off", widget=Dropdown(["off", "shift", "full"]))
    remove_internal_boundaries = Input(False, widget=CheckBox)
    non_destructive = Input(True, widget=CheckBox)
    use_obb = Input(False, widget=CheckBox)
    run_parallel = Input(True, widget=CheckBox)
    parallel_mode = Input(True, widget=CheckBox)
    fuzzy_value = Input(1.0e-7)
    check_inverted = Input(True, widget=CheckBox)

    to_keep: Optional[Sequence[BRep]] = Input(None)
    to_avoid: Optional[Sequence[BRep]] = Input(None)

    @Attribute
    def arguments(self):
        s1 = self.shape_in
        tools = self.tools
        if not s1:
            return tools
        else:
            lst = list(tools)
            lst.insert(0, s1)
            return lst

    def build(self):
        arguments = self.arguments

        if not isinstance(arguments, (list, tuple)):
            arguments = [arguments]

        # guard: will segfault otherwise
        if len(arguments) < 2:
            raise NotImplementedError("Not enough arguments to build something")

        this = BOPAlgo_CellsBuilder()
        # print("NonDestructive", this.NonDestructive())
        this.SetNonDestructive(self.non_destructive)
        # print("UseOBB", this.UseOBB())
        this.SetUseOBB(self.use_obb)
        # print("RunParallel", this.RunParallel())
        this.SetRunParallel(self.run_parallel)
        # print("GetParallelMode", this.GetParallelMode())
        this.SetParallelMode(self.parallel_mode)
        # print("FuzzyValue", this.FuzzyValue())
        this.SetFuzzyValue(self.fuzzy_value)
        # print("CheckInverted", this.CheckInverted())
        this.SetCheckInverted(self.check_inverted)

        pp_to_occ = {s: s.TopoDS_Shape for s in arguments}

        lst = create_toptools_listofshape(pp_to_occ.values())
        this.SetArguments(lst)
        glue = PY_TO_BOPAlgo_GlueEnum[self.glue]
        this.SetGlue(glue)
        this.Perform()

        to_keep = self.to_keep
        to_avoid = self.to_avoid
        if to_keep is None and to_avoid is None:
            this.AddAllToResult(1, False)
        else:
            if to_keep is None:
                to_keep = ()
            if to_avoid is None:
                to_avoid = ()
            to_keep = (pp_to_occ.pop(s) for s in to_keep)
            to_keep = create_toptools_listofshape(to_keep)
            to_avoid = (pp_to_occ.pop(s) for s in to_avoid)
            to_avoid = create_toptools_listofshape(to_avoid)

            this.AddToResult(to_keep, to_avoid, 1, False)

        if self.remove_internal_boundaries:
            this.RemoveInternalBoundaries()
        this.MakeContainers()
        return this


if __name__ == '__main__':
    from parapy.geom import Box, RectangularFace
    from parapy.gui import display
    arg1 = Box(1, 1, 1)
    arg2 = Box(1, 1, 1, position=arg1.position.translate(z=1))
    obj1 = GeneralFuse(tools=[arg1, arg2], glue="full")

    arg1 = RectangularFace(1, 1)
    arg2 = RectangularFace(1, 1,
                           position=arg1.position.translate(x=0.5, y=0.5))

    # mimics Fused(arg1, arg2)
    obj2 = GeneralFuse(tools=[arg1, arg2],
                       glue="shift")

    assert len(obj2.modified(arg1)) == 2
    assert len(obj2.modified(arg2)) == 2

    # mimics Subtracted(arg1, arg2)
    obj2 = GeneralFuse(tools=[arg1, arg2],
                       to_keep=[arg1],
                       to_avoid=[arg2],
                       glue="shift")

    assert len(obj2.modified(arg1)) == 1
    assert len(obj2.modified(arg2)) == 0

    # mimics Subtracted(arg2, arg1)
    obj2 = GeneralFuse(tools=[arg1, arg2],
                       to_keep=[arg2],
                       to_avoid=[arg1],
                       glue="shift")

    assert len(obj2.modified(arg1)) == 0
    assert len(obj2.modified(arg2)) == 1

    # Mimics Common(arg1, arg2)
    obj2 = GeneralFuse(tools=[arg1, arg2],
                       to_keep=[arg1, arg2],
                       glue="shift")

    assert len(obj2.modified(arg1)) == 1
    assert len(obj2.modified(arg2)) == 1

    # keep only arg1
    obj2 = GeneralFuse(tools=[arg1, arg2],
                       to_keep=[arg1],
                       glue="shift")

    assert len(obj2.modified(arg1)) == 2
    assert len(obj2.modified(arg2)) == 1

    # keep only arg2
    obj2 = GeneralFuse(tools=[arg1, arg2],
                       to_keep=[arg2],
                       glue="shift")

    assert len(obj2.modified(arg1)) == 1
    assert len(obj2.modified(arg2)) == 2

    arg1 = RectangularFace(2, 1)
    arg2 = RectangularFace(1, 1,
                           position=arg1.position.translate(x=1).rotate90('y'))
    obj3 = GeneralFuse(tools=[arg1, arg2])

    display((obj1, obj2, obj3))
