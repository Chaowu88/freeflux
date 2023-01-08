#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from sys import version_info
from warnings import warn
from setuptools import setup


if version_info[:2] < (3, 7):
    warn('FreeFlux is not tested in early versions of Python.')

setup()
