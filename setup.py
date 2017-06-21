from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy


ext_modules = [
    Extension('lunustbx',
                ['lunustbx.pyx'],
                extra_compile_args = ["-Ofast", "-fopenmp"],
                extra_link_args=['-fopenmp'],
                include_dirs=[numpy.get_include()])
]

setup(
    name = 'lunustbx Library',
    cmdclass = {'build_ext': build_ext},
    ext_modules=cythonize(ext_modules)
)