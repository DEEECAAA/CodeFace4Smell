#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='codeface',
    version='0.2.0',
    description='Codeface: Socio-Technical Analysis of Software Development',
    author='Wolfgang Mauerer',
    author_email='wolfgang.mauerer@oth-regensburg.de',
    url='https://github.com/siemens/codeface',
    packages=find_packages(exclude=['experiments']),
    package_data={'codeface': ['R/*.r', 'R/cluster/*.r', 'perl/*.pl']},
    entry_points={'console_scripts': ['codeface = codeface.cli:main']},
    install_requires=[
        'progressbar2',
        'pyctags',
        'PyYAML',
        'mysqlclient'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: OS Independent',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Quality Assurance'
    ],
    python_requires='>=3.6',
    license='GPLv2'
)
