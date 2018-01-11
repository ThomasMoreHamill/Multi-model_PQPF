from distutils.core import setup
from Cython.Build import cythonize
"""
python setup_calculate_quantile_for_precip_values.py build_ext --inplace
"""
setup(
    name = 'Calculate quantiles associated with precip amount for gamma distribution',
    ext_modules = cythonize("calculate_quantile_for_precip_values.pyx"),
)