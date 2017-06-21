# Diffuse_Scattering

This package takes raw diffraction images as input and returns 3D-reciprocal space maps of the diffuse X-ray scattering present in the images. The script relies on Lunus (https://github.com/mewall/lunus/tree/master), as well as DIALS methods within Phenix (https://www.phenix-online.org/download/nightly_builds.cgi).

## Quick Summary

### Contents

```
sematura_launcher.py        #script to run functions from library, can launch jobs on SGE
sematura.py                 #library of functions for processing diffraction images
sematura_params.py          #input parameters for diffuse data processing
lunustbx.pyx                #cython encoded functions from Lunus
setup.py                    #cythonize & compile *.pyx files
checkout.py                 #beginnings of a test file against gold standard outputs
*.npz                       #lattice files
*.pkl                       #crystal object from DIALS output
*.cbf                       #raw data to test pipeline
```

### Usage

```
libtbx.python setup.py build_ext --inplace      #one time only to compile lunustbx.pyx
libtbx.python sematura_launcher.py -i -p -a     #example of how to run pipeline in folder with raw data
```

## Detailed Instructions

```
Coming Soon!
```


### Prerequisities

```
DIALS
LUNUS
```


### Note

The script and the params file must be in the same directory for everything to run properly.

### Help

```
libtbx.python sematura_launcher.py -h
```




