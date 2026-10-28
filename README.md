# TLE SGP4 C++

A header-only C++17 TLE parser and near-Earth SGP4 propagator, written from Spacetrack Report #3 and Vallado et al. (2006).

It reproduces the reference SGP4 (the `sgp4` Python package, which wraps Vallado's code; WGS-72, improved mode) to under a micrometre in position over three days for the test cases in `tests/`, including the classic Spacetrack Report #3 TLE (88888, perigee 198 km, so it runs the simplified low-perigee drag branch) and an e = 0.18 orbit. Deep-space objects (period ≥ 225 min, which need SDP4) are detected and rejected instead of being propagated with the wrong model.

## Features

- TLE parsing (including the implied-decimal exponent fields) and epoch to Julian date
- Near-Earth SGP4: secular J2/J4 rates, atmospheric drag (BSTAR, C1-C5, D2-D4), long- and short-period periodics
- State vectors in TEME [km, km/s]; GMST for rotating to Earth-fixed
- Period, semi-major axis, apogee/perigee altitudes
- Single header, no dependencies

## Quick Start

```cpp
#include "sgp4_propagator.hpp"

using namespace sgp4;

int main() {
    // ISS TLE
    std::string line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025";
    std::string line2 = "2 25544  51.6400 208.9163 0006703 296.7361  63.2861 15.49560921426457";

    // Initialize SGP4
    SGP4 sgp4(line1, line2, "ISS");

    // Propagate 30 minutes from epoch
    StateVector sv = sgp4.propagate(30.0);

    std::cout << "Position: " << sv.x << ", " << sv.y << ", " << sv.z << " km\n";
    std::cout << "Velocity: " << sv.vx << ", " << sv.vy << ", " << sv.vz << " km/s\n";

    return 0;
}
```

## Building

### Requirements
- C++17 compatible compiler
- CMake 3.14+ (optional)

### With CMake

```bash
mkdir build && cd build
cmake ..
make

# Run examples and the reference test
./iss_tracking
./constellation_demo
ctest --output-on-failure
```

### Direct Compilation

```bash
g++ -std=c++17 -O2 -I include examples/iss_tracking.cpp -o iss_tracking
```

## Examples

### ISS Tracking (`examples/iss_tracking.cpp`)

Tracks the International Space Station and shows its ground track over one orbit. Displays orbital elements, computes position/velocity states, and shows how the altitude varies.

### Constellation Demo (`examples/constellation_demo.cpp`)

Loads several TLEs (ISS, Hubble, Starlink, plus a GPS and a Molniya satellite that are skipped as deep-space), compares their orbits, tracks altitude over 24 hours, and finds the closest approach to the ISS.

### Plotting (`scripts/plot_examples.py`)

Python script that generates visualizations from the example outputs. Run `python3 scripts/plot_examples.py` after building to create plots.

## API Reference

### TLE Structure

```cpp
struct TLE {
    std::string name;
    int satellite_number;
    char classification;
    int epoch_year;
    double epoch_day;
    double mean_motion_dot;
    double mean_motion_ddot;
    double bstar;
    double inclination;      // [deg]
    double raan;             // [deg]
    double eccentricity;
    double arg_perigee;      // [deg]
    double mean_anomaly;     // [deg]
    double mean_motion;      // [rev/day]
    double epoch_jd;         // Julian Date
};
```

### StateVector Structure

```cpp
struct StateVector {
    double x, y, z;      // Position [km] in TEME
    double vx, vy, vz;   // Velocity [km/s] in TEME
    
    double position_magnitude() const;
    double velocity_magnitude() const;
};
```

### SGP4 Class

```cpp
class SGP4 {
public:
    // Constructors
    SGP4(const TLE& tle);
    SGP4(const std::string& line1, const std::string& line2, 
         const std::string& name = "");

    // Propagation
    StateVector propagate(double tsince);    // minutes since epoch
    StateVector propagateJD(double jd);      // to Julian Date

    // Orbital parameters
    double getSemiMajorAxis() const;  // [km]
    double getPeriod() const;         // [min]
    double getApogee() const;         // [km altitude]
    double getPerigee() const;        // [km altitude]
    bool isDeepSpace() const;
};
```

### Utility Functions

```cpp
// Parse TLE from strings
TLE parseTLE(const std::string& line1, const std::string& line2, 
             const std::string& name = "");

// Time conversions
double epochToJD(int year, double day);
std::string jdToDateString(double jd);

// Output
void printStateVector(const StateVector& sv, const std::string& label = "");
```

## TLE Format

Two-Line Element sets are NORAD's standard format for sharing orbital data:

```
ISS (ZARYA)
1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025
2 25544  51.6400 208.9163 0006703 296.7361  63.2861 15.49560921426457
```

**Line 1** has the satellite ID, epoch time, drag terms, and BSTAR coefficient.  
**Line 2** contains the Keplerian elements: inclination, RAAN, eccentricity, argument of perigee, mean anomaly, and mean motion.

## Coordinate Frame

Outputs are in the **TEME (True Equator Mean Equinox)** frame, a quasi-inertial Earth-centered frame. To get latitude/longitude for ground tracks, you'll need to rotate by Greenwich Sidereal Time (GMST). Check `examples/iss_tracking.cpp` for a working example.

## Validation

`tests/test_sgp4.cpp` compares position and velocity with the reference implementation at 0, 6, 24 and 72 hours for three TLEs, and also checks the epoch Julian date and the deep-space rejection:

```
max position difference vs reference: 6.6e-07 m
max velocity difference vs reference: 6.5e-07 mm/s
```

`scripts/make_reference.py` regenerates the reference states (`pip install sgp4`).

Agreement with the reference only shows that the code implements SGP4 correctly. SGP4 itself is a mean-element model fit to tracking data: for LEO it is typically good to around a kilometre near the TLE epoch, and the error grows by kilometres per day, so use recent TLEs.

## Limitations

- **Near-Earth only.** SDP4 (lunar/solar perturbations and the 12-h/24-h resonances for GPS, Molniya, GEO) is not implemented. `propagate()` throws `std::domain_error` for those objects.
- Output is TEME. GMST is the IAU-82 formula with UT1 ≈ UTC. No polar motion or TEME→GCRF conversion is included.
- Errors that the reference code flags (decay, eccentricity out of range) are raised as `std::domain_error`.

## Where to Get TLEs

- [Space-Track.org](https://www.space-track.org) (free account)
- [CelesTrak](https://celestrak.org) (no registration)

## References

1. F. R. Hoots and R. L. Roehrich, *Spacetrack Report No. 3: Models for Propagation of NORAD Element Sets*, 1980.
2. D. A. Vallado, P. Crawford, R. Hujsak and T. S. Kelso, "Revisiting Spacetrack Report #3," AIAA 2006-6753.
3. B. Rhodes, `sgp4` Python package (reference used for the tests).

## License

MIT. See [LICENSE](LICENSE).


