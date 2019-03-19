.. lightcurve documentation master file, created by
   sphinx-quickstart on Tue Mar 19 11:09:30 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to lightcurve's documentation!
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


**lightcurve** is a module that creates light curves of a similar field of stars.
This is useful in time domain astronomy in search of transient events
and studies of variable stars, i.e. eclipsing binaries.

Usage
-----
Simple to use:

.. code-block:: python

 import lightcurves as lc

 object_table = lc.makelist('')

 # parallelized alignment
 object_table = lc.paralign(object_table)

 # or serial alignment
 object_table = lc.align(object_table)

 # All-in-one
 object_table = lc.do_lightcurve(object_table)

Process
-------
The process follows:

1. Take image list

 a. Options (or considerations)

  i. Reduce images

  ii. Assign WCS

2. Align images (How does reducing affect WCS?)

 a. [Create a deep sky image pending]

3. Source extraction [from deepsky image]

4. Match sources between images [pending]

5. Create Timeseries **plot** of different sources

6. Create **reference image** with all candidate sources circled

Contribute
----------
- Issue Tracker: https://github.com/moemyself3/lightcurve/issues
- Source Code: https://github.com/moemyself3/lightcurve

License
-------
This project is licensed under the MIT license.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
