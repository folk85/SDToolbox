from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(name='SDToolbox',
    version='0.0.1',
    description='The Shock & Detonation Toolbox',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/folk85/SDToolbox',
    project_urls={
        'Documentation': 'http://shepherd.caltech.edu/EDL/public/cantera/html/SD_Toolbox/',
        # 'Funding': 'https://donate.pypi.org',
        # 'Say Thanks!': 'http://saythanks.io/to/example',
        'Source': 'https://github.com/folk85/SDToolbox/',
        'Tracker': 'https://github.com/folk85/SDToolbox/issues',
    },
    author='Sergey Medvedev',
    author_email='medvsn@gmail.com',
    license='MIT',
    packages=['SDToolbox'],
    python_requires='>2.6, <4',
    install_requires=[
          'scipy',
          'numpy',
          'cantera',
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
    include_package_data=True,
    zip_safe=False)