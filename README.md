# Diffuse_Scattering

This script takes raw diffraction images as input and returns 3D-reciprocal space maps of the diffuse X-ray scattering present in the images. The script relies on the Lunus (https://github.com/mewall/lunus/tree/master), as well as DIALS methods within Phenix (https://www.phenix-online.org/download/nightly_builds.cgi).

## Contents

sematura.py                 wrapper for Lunus & DIALS software for processing diffraction images
sematura_params.py          input parameters for diffuse data processing

### Prerequisities

You need to have Phenix installed and properly sourced.


You need to download and compile the Lunus software package.

* Download Lunus

```
$ git clone -b master https://github.com/mewall/lunus.git
```

* Move Lunus, if desired

```
$ mv ../path_to/Downloads/lunus ../new_path_to/lunus
```

* Note for Mac Users

Make sure you have gcc installed, currently in osx the gcc command actually runs clang, which does not support openMP for parallel processing (i.e.:lunus will not compile properly)

* Compile Lunus

```
$ cd ../path_to/lunus/c/src

$ pwd

$ cat 00README

$ csh

% setenv C_HOME path_to/lunus/c

% makemake Makefile
```
Mac Users: Open Makefile with text editor of your choice, and change “gcc” on line 3 to your version (eg: “gcc-6”) then save the file

```
% make all
```

You can exit csh now.


### Installing

To run sematura.py you may need to make the following changes to Phenix.

Open the following script in the text editor of your choice:
/path_to/phenix_version/modules/cctbx_project/iotbx/detectors/detectorbase.py

* add “+” lines and comment out “-“ lines

```
detectorbase.py
===================================================================
--- detectorbase.py (revision 18023)
+++ detectorbase.py (working copy)
@@ -203,9 +203,7 @@
    F.close()
    from iotbx.detectors import WriteADSC
    if mod_data==None: mod_data=self.linearintdata
-    if not mod_data.all_ge(0):
-      from libtbx.utils import Sorry
-      raise Sorry("Negative values not allowed when writing SMV")
+    mod_data = mod_data.set_selected(mod_data < 0, 0)
    WriteADSC(fileout,mod_data,self.size1,self.size2,endian)
```

