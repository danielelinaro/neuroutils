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

def _test(verbosity=2):

    from .test_node import TestNode
    from .test_tree import TestTree
    from unittest import TestSuite, TestLoader, TextTestRunner

    # list of test cases
    suites = [
        TestLoader().loadTestsFromTestCase(TestNode),
        TestLoader().loadTestsFromTestCase(TestTree),
    ]
    TextTestRunner(verbosity=verbosity).run(TestSuite(suites))
