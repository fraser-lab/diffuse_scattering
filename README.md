# Diffuse_Scattering

This package takes raw diffraction images as input and returns 3D-reciprocal space maps of the diffuse X-ray scattering present in the images. The script relies on a custom build of DIALS (http://dials.lbl.gov/installation.html), which incorporates methods from the Lunus software package (https://github.com/mewall/lunus/tree/master).

## Quick Summary

### Contents

```
sematura_launcher.py        #script to run functions from library, can launch jobs on SGE
sematura.py                 #library of functions for processing diffraction images
sematura_params.py          #input parameters for diffuse data processing
checkout.py                 #beginnings of a test file against gold standard outputs
*.npz                       #lattice files
*.pkl                       #crystal object from DIALS output
*.cbf                       #raw data to test pipeline
```

### Usage

```
libtbx.python sematura_launcher.py -i -p -a     #example of how to run pipeline in folder with raw data
```

## Detailed Instructions

```
Coming Soon!
```


### Prerequisities

```
DIALS (Custom Build)
```


### Help

```
libtbx.python sematura_launcher.py -h
```




