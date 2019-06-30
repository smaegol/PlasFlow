#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools


setuptools.setup(
    name='plasflow',
    version='1.2.0',
    description="Identify plasmids using a nueral net. Original work by Pawel Krawczyk refactored by David Danko",
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
        'biopython',
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
