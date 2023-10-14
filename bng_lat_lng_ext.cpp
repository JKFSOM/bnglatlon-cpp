#define _USE_MATH_DEFINES

#include <Python.h>
#include <cmath>
#include <tuple>

inline double secondsToRadians(double sec)
{
    return sec * M_PI / (180.0 * 3600.0);
}

std::tuple<double, double> OSGB36toWGS84(double E, double N)
{
    // The Airy 1830 semi-major and semi-minor axes used for OSGB36 (m)
    constexpr double a = 6377563.396;
    constexpr double b = 6356256.909;

    // scale factor on the central meridian
    constexpr double F0 = 0.9996012717;

    // Latitude and longtitude of true origin (radians)
    double lat0 = 49.0 * M_PI / 180.0;
    double lon0 = -2.0 * M_PI / 180.0;

    // Northing & easting of true origin (m)
    constexpr double N0 = -100000.0;
    constexpr double E0 = 400000.0;

    // eccentricity squared
    const double e2 = 1.0 - (b * b) / (a * a);

    const double n = (a - b) / (a + b);

    // Initilaise the iterative variables
    double lat = lat0;
    double M = 0.0;

    while (N - N0 - M >= 0.00001)
    {
        lat = (N - N0 - M) / (a * F0) + lat;

        double M1 = (1.0 + n + (5.0 / 4.0) * pow(n, 2.0) + (5.0 / 4.0) * pow(n, 3.0)) * (lat - lat0);
        double M2 = (3.0 * n + 3.0 * pow(n, 2.0) + (21.0 / 8.0) * pow(n, 3.0)) * sin(lat - lat0) * cos(lat + lat0);
        double M3 = ((15.0 / 8.0) * pow(n, 2.0) + (15.0 / 8.0) * pow(n, 3.0)) * sin(2.0 * (lat - lat0)) * cos(2.0 * (lat + lat0));
        double M4 = (35.0 / 24.0) * pow(n, 3.0) * sin(3.0 * (lat - lat0)) * cos(3.0 * (lat + lat0));

        // meridional arc
        M = b * F0 * (M1 - M2 + M3 - M4);
    }

    const double sin_lat_squared = pow(sin(lat), 2.0);
    // transverse radius of curvature
    double nu = a * F0 / sqrt(1.0 - e2 * sin_lat_squared);

    // meridional radius of curvature
    double rho = a * F0 * (1.0 - e2) * pow(1.0 - e2 * sin_lat_squared, -1.5);
    double eta2 = nu / rho - 1.0;

    double sec_lat = 1.0 / cos(lat);

    double VII = tan(lat) / (2.0 * rho * nu);
    double VIII = tan(lat) / (24.0 * rho * pow(nu, 3.0)) * (5.0 + 3.0 * pow(tan(lat), 2.0) + eta2 - 9.0 * pow(tan(lat), 2.0) * eta2);
    double IX = tan(lat) / (720.0 * rho * pow(nu, 5.0)) * (61.0 + 90.0 * pow(tan(lat), 2.0) + 45.0 * pow(tan(lat), 4.0));
    double X = sec_lat / nu;
    double XI = sec_lat / (6.0 * pow(nu, 3.0)) * (nu / rho + 2.0 * pow(tan(lat), 2.0));
    double XII = sec_lat / (120.0 * pow(nu, 5.0)) * (5.0 + 28.0 * pow(tan(lat), 2.0) + 24.0 * pow(tan(lat), 4.0));
    double XIIA = sec_lat / (5040.0 * pow(nu, 7.0)) * (61.0 + 662.0 * pow(tan(lat), 2.0) + 1320.0 * pow(tan(lat), 4.0) + 720.0 * pow(tan(lat), 6.0));

    double dE = (E - E0);
    // These are on the wrong ellipsoid currently: Airy1830 (denoted by _1)
    double lat1 = lat - VII * pow(dE, 2.0) + VIII * pow(dE, 4.0) - IX * pow(dE, 6.0);
    double lon1 = lon0 + X * dE - XI * pow(dE, 3.0) + XII * pow(dE, 5.0) - XIIA * pow(dE, 7.0);

    // Convert to GRS80 ellipsoid
    // First convert to cartesian from spherical polar coordinates
    double H = 0.0; // Third spherical coord.
    double x1 = (nu / F0 + H) * cos(lat1) * cos(lon1);
    double y1 = (nu / F0 + H) * cos(lat1) * sin(lon1);
    double z1 = ((1.0 - e2) * nu / F0 + H) * sin(lat1);

    // Perform Helmert transform (to go between Airy 1830 and GRS80)
    double s = -20.4894 * pow(10.0, -6.0); // The scale factor -1

    // Translation along x,y, and z axes respectively
    constexpr double tx = 446.448;
    constexpr double ty = -125.157;
    constexpr double tz = 542.060;

    // Rotate About x-axis, y-axis, and z-axis respectively
    constexpr double rxs = 0.1502;
    constexpr double rys = 0.2470;
    constexpr double rzs = 0.8421;

    // convert seconds to radians
    const double rx = secondsToRadians(rxs);
    const double ry = secondsToRadians(rys);
    const double rz = secondsToRadians(rzs);

    // perform transformation
    const double x2 = tx + (1.0 + s) * x1 + (-rz * y1) + (ry * z1);
    const double y2 = ty + (rz * x1) + (1.0 + s) * y1 + (-rx * z1);
    const double z2 = tz + (-ry * x1) + (rx * y1) + (1.0 + s) * z1;

    // Back to spherical polar coordinates from cartesian
    // Need some of the characteristics of the new ellipsoid

    // The GSR80 semi-major and semi-minor axes used for WGS84(m)
    constexpr double a2 = 6378137.000;
    constexpr double b2 = 6356752.3141;

    // the eccentricity of the GRS80 ellipsoid
    double e2_2 = 1.0 - (b2 * b2) / (a2 * a2);
    double p = sqrt(pow(x2, 2.0) + pow(y2, 2.0));

    // Lat is obtained by an iterative proceedure:
    lat = atan2(z2, (p * (1.0 - e2_2))); // Initial value
    double latold = 2.0 * M_PI;

    double nu2;
    while (abs(lat - latold) > pow(10.0, -16.0))
    {
        // Make a copy of lat, because we're about to overwrite it
        double lat_copy = lat;
        lat = latold;
        latold = lat_copy;

        nu2 = a2 / sqrt(1.0 - e2_2 * pow(sin(latold), 2.0));
        lat = atan2(z2 + e2_2 * nu2 * sin(latold), p);
    }

    // Lon and height are then pretty easy
    double lon = atan2(y2, x2);
    H = p / cos(lat) - nu2;

    // Convert to degrees
    lat = lat * 180.0 / M_PI;
    lon = lon * 180.0 / M_PI;

    double lat_rounded = round(lat * 1000000.0) / 1000000.0;
    double lon_rounded = round(lon * 1000000.0) / 1000000.0;

    return std::make_tuple(lat_rounded, lon_rounded);
}

static PyObject *py_OSGB36toWGS84(PyObject *self, PyObject *args)
{

    double E, N;

    if (!PyArg_ParseTuple(args, "dd", &E, &N))
        return NULL;

    double lat, lon;
    std::tie(lat, lon) = OSGB36toWGS84(E, N);

    return Py_BuildValue("dd", lat, lon);
}

static PyMethodDef OSGB36toWGS84Methods[] = {
    {"OSGB36toWGS84", py_OSGB36toWGS84, METH_VARARGS, "Convert OSGB36 coordinates to WGS84."},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef OSGB36toWGS84_module = {
    PyModuleDef_HEAD_INIT,
    "OSGB36toWGS84_ext",
    NULL,
    -1,
    OSGB36toWGS84Methods};

PyMODINIT_FUNC PyInit_OSGB36toWGS84_ext(void)
{

    return PyModule_Create(&OSGB36toWGS84_module);
}