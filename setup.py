from setuptools import setup

# Get the version from file itself (not imported)
with open('lightcurve.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            _, _, lc_version = line.replace("'", '').split()
            break

with open('README.md', 'r') as f:
    long_description = f.read()

setup(name='lightcurve',
      version=lc_version,
      description='Creates lightcurve database of series of images.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Moises Castillo',
      author_email='castillo.moises11@gmail.com',
      url='https://github.com/moemyself3/lightcurve',
      py_modules=['lightcurve', ],
      install_requires=["astroalign==1.0.1",
                        "astropy==3.1.2",
                        "matplotlib==2.0.2",
                        "numpy==1.14.2",
                        "photutils==0.6",
                        ],
      test_suite='tests',
      )
