#!/usr/bin/env python
""" IMaRS WV2 processing setup.py """

from setuptools import setup
import io

VERSION = '0.0.1'  # should match __version__ in PROJECT_NAME.__init__.py


# === long_description comes from README.md
def read(*filenames, **kwargs):
    """Helper fn to help get project info for setup()"""
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')  # , 'CHANGES.txt')

# === installation requirements read from requirements files
# so these two commands are equivalent:
# * `pip install .`
# * `pip install -r requirements.txt`
_tests_require = [
    line.strip() for line in open('requirements_tests.txt')
    if line.strip() and not line.strip().startswith('--')
]

_install_requires = [
    line.strip() for line in open('requirements.txt')
    if line.strip() and not line.strip().startswith('--')
]

_extras_require = {  # allows install w/ `pip install .[test]`
    'test': _tests_require
}

setup(
    name='wv_classify',
    description='IMaRS WV2 classification processing scripts',
    author='Matthew McCarthy',
    author_email='mjm8@mail.usf.edu',
    url='https://github.com/usf-imars/wv2-processing',
    version=VERSION,
    long_description=long_description,
    install_requires=_install_requires,
    tests_require=_tests_require,
    extras_require=_extras_require,
)
