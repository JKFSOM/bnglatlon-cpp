from setuptools import setup, Extension

bng_lat_lng_module = Extension(
    "OSGB36toWGS84_ext",
    sources=["bng_lat_lng_ext.cpp"],
    extra_compile_args=["-std=c++17"],
)

setup(
    name="BngLatLngExtension",
    version="0.1",
    description="Python extension in C++ to convert BNG to LatLng",
    ext_modules=[bng_lat_lng_module],
)
