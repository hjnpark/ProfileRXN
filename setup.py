from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='ProfileRXN',
    description='Given two xyz files, it will interpolate and calculate energies of reactant, product and transition state.',
    url='https://github.com/hjnpark/ProfileRXN',
    author='Heejune Park',
    packages=find_packages(),
    package_data={'': ['*.ini']},
    include_package_data=True,
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={'console_scripts': [
        'profile-rxn = profilerxn.profilerxn:main',
    ]},
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
    zip_safe=True,
)
