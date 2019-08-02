[![Documentation Status](https://readthedocs.org/projects/lightcurator/badge/?version=latest)](https://lightcurator.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/124146048.svg)](https://zenodo.org/badge/latestdoi/124146048)


# LightCurve

LightCurve is a module that creates light curves of a similar field of stars.
This is useful in time domain astronomy in search of transient events
and studies of variable stars, i.e. eclipsing binaries.


# Table of Contents

# Installation

# Usage

```
from lightcurator import lightcurves as lc

object_table = lc.makelist('')

# parallelized alignment
object_table = lc.paralign(object_table)

# or serial alignment
object_table = lc.align(object_table)

# All-in-one
object_table = lc.do_lightcurve(object_table)
```

The process follows:
1. Take image list
2. Align images
3. Create a deepsky image
4. Plate solve deepsky using `astrometery`
5. Sorce extraction from deepsky image
6. Create **reference image** with all candidate sources circled
7. Match sources between aligned images
8. Cross match sources from aligned images with sources from deepsky image
9. Create Timeseries **plot** of different sources
10. Cross match sources with catalogs like VSX and GCVS 


# Credits

# License

