#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2021 ParaPy Holding B.V.
#
# You may use the contents of this file in your application code.
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR
# PURPOSE.

# from parapy.lib.code_aster.primitive.material import Material


class Material:
    def __init__(self, name, E, nu, density):
        self.name = name
        self.E = E
        self.nu = nu
        self.density = density


STEEL = Material(name="steel", E=3e11, nu=0.2, density=7850)
ALUMINIUM = Material(name="aluminium", E=7e10, nu=0.33, density=2700)

