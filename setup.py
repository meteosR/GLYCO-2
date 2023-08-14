from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules = [
    Extension(
        r'bundle_cylinder_calc',
        [r'bundle_cylinder_calc.pyx']
    ),
]

setup(
    name='bundle_cylinder_calc',
    ext_modules=cythonize(ext_modules),
)

# python3 setup.py build_ext --inplace