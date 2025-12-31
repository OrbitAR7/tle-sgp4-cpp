# TLE SGP4 C++

A simple header-only C++ library for satellite orbit propagation using SGP4. Perfect for tracking satellites, predicting passes, or analyzing constellations. No external dependencies needed - just drop the header in and you're good to go.

## Features

- **TLE Parsing** - Handles standard two-line element sets
- **SGP4 Propagation** - Near-Earth orbit propagation (includes J2-J4 perturbations)
- **Position & Velocity** - Returns state vectors in TEME frame
- **Orbital Parameters** - Calculates period, apogee/perigee, semi-major axis
- **Time Conversions** - Julian Date support built-in
- **Header-only** - Just one file, no dependencies or linking needed

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

# Run examples
./iss_tracking
./constellation_demo

# Run tests
./test_sgp4
```

### Direct Compilation

```bash
g++ -std=c++17 -O2 -I include examples/iss_tracking.cpp -o iss_tracking
```

## Examples

### ISS Tracking (`examples/iss_tracking.cpp`)

Tracks the International Space Station and shows its ground track over one orbit. Displays orbital elements, computes position/velocity states, and shows how the altitude varies.

### Constellation Demo (`examples/constellation_demo.cpp`)

Tracks multiple satellites with different orbit types (LEO, MEO, HEO) over 24 hours. Really interesting to see how a Molniya orbit's altitude swings compared to circular orbits like the ISS or GPS.

<img width="1400" height="800" alt="constellation_altitude" src="https://github.com/user-attachments/assets/d2ef93a4-1a84-442d-ba6a-edb251a663ec" />

*24-hour altitude evolution showing different orbit types - notice how the Molniya (red) swings from ~1,600 km up to ~38,000 km*

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

Outputs are in the **TEME (True Equator Mean Equinox)** frame - basically an Earth-centered inertial frame. To get latitude/longitude for ground tracks, you'll need to rotate by Greenwich Sidereal Time (GMST). Check `examples/iss_tracking.cpp` for a working example.

## Limitations

This is designed for **near-Earth satellites** (orbital period < 225 minutes). Deep-space objects like GEO or Molniya use a simplified propagation model here - it works but isn't as accurate as a full SDP4 implementation.

If you need production-grade deep-space propagation, check out Vallado's reference code or the official Space-Track libraries.

## Accuracy

SGP4 is pretty good but not perfect:
- **Position:** 1-3 km error at epoch
- **Velocity:** 1-3 m/s error at epoch  
- Degrades about 1-3 km per day from the epoch

For best results, grab fresh TLEs (< 2 days old) from Space-Track or Celestrak.

## Where to Get TLEs

- [Space-Track.org](https://www.space-track.org) - Need to register but it's free
- [Celestrak](https://celestrak.org) - No registration required  
- [N2YO](https://www.n2yo.com) - Easy to browse and search

## References

Based on the classic SGP4 papers and Vallado's revisit:

1. Hoots & Roehrich (1980) - *Spacetrack Report No. 3* (the original SGP4 spec)  
2. Vallado et al. (2006) - *Revisiting Spacetrack Report #3* (modern updates)  
3. [Celestrak SGP4 Resources](https://celestrak.org/software/tskelso-sw.php)

## License

MIT License — See [LICENSE](LICENSE) for details.


