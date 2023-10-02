from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        r'bundle_cylinder_calc',
        [r'bundle_cylinder_calc.pyx'],
    ),
]

setup(
    name='bundle_cylinder_calc',
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()]
)

# python3 setup.py build_ext --inplace