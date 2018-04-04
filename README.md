# LightCurve

LightCurve is a module that creates light curves of a similar field of stars.
This is useful in time domain astronomy in search of transient events
and studies of variable stars, i.e. eclipsing binaries.


# Table of Contents

# Installation

# Usage

```
import lightcurves as lc

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
   - Options (*or considerations)
     - Reduce images
     - Assign WCS*
2. Align images (How does reducing affect WCS?)
2.1 [Create a deep sky image pending]
3. Sorce extraction [from deepsky image]
4. Match sources between images [pending]
5. Create Timeseries **plot** of different sources
6. Create **reference image** with all candidate sources circled

# Credits

# License

