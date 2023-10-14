import OSGB36toWGS84_ext as OSGB36toWGS84_cpp

from bng_latlon import OSGB36toWGS84
import timeit

EASTING = 51.4778
NORTHING = -0.0014


def cpp():
    """
    C++ version of the function
    """
    _, _ = OSGB36toWGS84_cpp.OSGB36toWGS84(EASTING, NORTHING)


def py():
    """
    Python version of the function
    """
    _, _ = OSGB36toWGS84(EASTING, NORTHING)


if __name__ == "__main__":
    RUNS = 100

    cpp_time = timeit.timeit(cpp, number=RUNS)
    print("cpp:", cpp_time, "seconds")
    python_time = timeit.timeit(py, number=RUNS)
    print("python:", python_time, "seconds")
    delta = python_time - cpp_time
    print("delta:", delta, "seconds")
    delta_percentage = (cpp_time - python_time) / python_time * 100
    print("delta percentage (%): ", delta_percentage)
