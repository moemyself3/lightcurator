from setuptools import setup, find_packages

# Get the version from file itself (not imported)
with open('lightcurator/__init__.py', 'r') as f:
    for line in f:
        if line.startswith("__version__"):
            _, _, lc_version = line.replace("'", '').split()


with open('README.md', 'r') as f:
    long_description = f.read()

setup(name='lightcurator',
    version=lc_version,
    description='Creates a lightcurve database of series of images.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Moises Castillo',
    author_email='castillo.moises11@gmail.com',
    url='https://github.com/moemyself3/lightcurator',
    packages=find_packages(),
    install_requires=['astroalign==1.0.1',
                      'astropy>=3.1.2',
                      'astroquery>=0.3',
                      'matplotlib>=2.0.2',
                      'numpy>=1.14.2',
                      'photutils>=0.6',
                      'ccdproc>=1.2.0'
                      ],
    test_suite='tests',
    license='MIT License',
    classifiers=[
        'License :: OSI Approved :: MIT License'
        ]
     )
