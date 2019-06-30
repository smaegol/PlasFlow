#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools


setuptools.setup(
    name='plasflow',
    version='1.1.0',
    description="",
    author="Pawel Krawczyk",
    author_email='p.krawczyk@ibb.waw.pl',
    url='https://github.com/smaegol/PlasFlow',
    packages=setuptools.find_packages(),
    package_dir={'plasflow': 'plasflow'},
    install_requires=[
        'argparse',
        'pandas',
        'scipy',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'plasflow=plasflow.cli:cli_main'
        ]
    },
    classifiers=[
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3',
    ],
    package_data={'models': ['models/*']},
)
