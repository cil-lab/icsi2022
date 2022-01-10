from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


icsi2022_extension = Extension(
    name="icsi2022",
    sources=["icsi2022.pyx"],
    libraries=["icsi2022"],
    library_dirs=["lib"],
    include_dirs=["lib"],
)
setup(name="icsi2022", ext_modules=cythonize([icsi2022_extension]), include_dirs=[numpy.get_include()], compiler_directives={'language_level' : "3"})
