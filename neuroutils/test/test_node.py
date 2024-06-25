# -*- coding: utf-8 -*-
"""
Copyright (C) 2024 Daniele Linaro <danielelinaro@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""


import unittest
from ..tree import Node

class TestNode(unittest.TestCase):

    def test_equal(self):
        A,B = Node(ID='0'), Node(ID='0')
        self.assertEqual(A, B, 'Nodes are not equal')

    def test_not_equal(self):
        A,B = Node(ID='0'), Node(ID='1')
        self.assertNotEqual(A, B, 'Nodes are equal')

