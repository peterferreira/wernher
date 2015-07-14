from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import os
import sys

import wernher

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md') #, 'CHANGES.md')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='wernher',
    version=wernher.__version__,
    url='https://github.com/theodoregoetz/kepler-kerman',
    license='GNU GPL V3',
    author='Johann Goetz',
    tests_require=['pytest'],
    install_requires=['numpy>=1.8.0',
                    'scipy>=0.14.1',
                    'krpc>=0.1.9',
                    ],
    cmdclass={'test': PyTest},
    author_email='theodore.goetz@gmail.com',
    description='Toolkit for Kerbal Space Program with kRPC',
    long_description=long_description,
    packages=['wernher'],
    include_package_data=True,
    platforms='any',
    test_suite='wernher.test',
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
        'Natural Language :: English',
        'Environment :: Console',
        'Environment :: X11 Applications',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Games/Entertainment :: Simulation',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
    extras_require={
        'testing': ['pytest'],
    }
)
