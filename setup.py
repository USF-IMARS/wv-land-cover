#!/usr/bin/env python
""" IMaRS WV2 processing setup.py """

from setuptools import setup
import io

VERSION='0.0.1'

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md') #, 'CHANGES.txt')

setup(name='wv2_processing',
    version=VERSION,
    description='IMaRS WV2 classification processing scripts',
    long_description=long_description,
    author='Matthew McCarthy',
    author_email='mjm8@mail.usf.edu',
    url='https://github.com/usf-imars/wv2-processing',
    tests_require=['nose'],
    install_requires=['gdal']
)
