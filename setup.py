import os
import re
from setuptools import find_packages

try:
    from setuptools import setup
    from distutils.core import setup
except ImportError:
    pass

version = __import__('svkits').get_version()

_REQUIREMENTS_FILE = 'REQUIREMENTS.txt'
_README = 'README.md'


def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


def _get_description(file_name):
    with open(file_name, 'r') as f:
        _long_description = f.read()
    return _long_description


def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    rx = re.compile('^[A-z]')
    requirements = [l for l in lines if rx.match(l) is not None]
    if "READTHEDOCS" in os.environ:
        requirements = [r for r in requirements if not "pbcore" in r]
    return requirements

setup(
    name='svkits',
    version=version,
    package_dir={'svkits': 'svkits'},
    packages=find_packages('.'),
    license='BSD',
    author='yli',
    author_email='yli@pacificbiosciences.com',
    description='Tool kits for structure variation simulation, validation, analysis, tests and ect.',
    setup_requires=['nose>=1.0'],
    # Maybe the pbtools-* should really be done in a subparser style
    entry_points={'console_scripts': [
        'make-subreads-bam-from-fasta = svkits.make_subreads_bam_from_fasta:main',
        'make-subreads-bam-from-zmws = svkits.make_subreads_bam_from_zmws:main',
        'validate-sv-bed = svkits.validate_sv_bed:main',
        'add-an-indel-to-fasta = svkits.add_an_indel_to_fasta:main',
        'add-indels-to-fasta = svkits.add_indels_to_fasta:main',
        'mummer-plot = svkits.mummer_plot:main',
        'sim-x2coverage = svkits.sim_x2coverage:main'
    ]},
    install_requires=_get_requirements(_get_local_file(_REQUIREMENTS_FILE)),
    tests_require=['nose'],
    long_description=_get_description(_get_local_file(_README)),
    classifiers=['Development Status :: 4 - Beta'],
    include_package_data=True,
    zip_safe=False
)
