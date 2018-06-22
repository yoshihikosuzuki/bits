# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages
import sys


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('requirements.txt') as f:
    requirements = f.read().strip().split('\n')
    if sys.version_info[0] == 3:
        print('\033[92m' + '[INFO] ]Install pbcore from https://github.com/yoshihikosuzuki/pbcore if you are using python3.' + '\033[0m')

setup(
    name='BITS',
    version='0.0.1',
    description='BioInformatics ToolS',
    long_description=readme,
    author='Yoshihiko Suzuki',
    author_email='ys.neoteny@gmail.com',
    install_requires=requirements,
    url='https://github.com/yoshihikosuzuki/BITS',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
