/**
 * ISS Tracking Example
 * 
 * This example shows how to track the ISS using SGP4 propagation.
 * We parse the TLE, compute ground tracks, and visualize the orbit over time.
 * 
 * Author: Hamoud Alshammari
 * 2025
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "../include/sgp4_propagator.hpp"

using namespace sgp4;

// Convert TEME coordinates to lat/lon/alt for ground track visualization.
// Uses GMST rotation (ignores polar motion, which is fine for plotting).
struct GroundTrack {
    double lat;   // degrees
    double lon;   // degrees  
    double alt;   // km
};

GroundTrack temeToGroundTrack(const StateVector& sv, double jd) {
    using namespace constants;

    // TEME -> pseudo-ECEF via GMST
    double gst = SGP4::gstime(jd);
    double cos_gst = std::cos(gst);
    double sin_gst = std::sin(gst);
    double x_ecef = sv.x * cos_gst + sv.y * sin_gst;
    double y_ecef = -sv.x * sin_gst + sv.y * cos_gst;
    double z_ecef = sv.z;

    // WGS-84 ellipsoid
    constexpr double a = 6378.137;                 // km
    constexpr double f = 1.0 / 298.257223563;
    constexpr double b = a * (1.0 - f);
    constexpr double e2 = 1.0 - (b * b) / (a * a);

    double lon = std::atan2(y_ecef, x_ecef);
    double r_xy = std::hypot(x_ecef, y_ecef);
    double lat = std::atan2(z_ecef, r_xy); // initial guess

    // Iterate latitude for ellipsoid height
    for (int i = 0; i < 5; ++i) {
        double sin_lat = std::sin(lat);
        double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
        double alt = r_xy / std::cos(lat) - N;
        lat = std::atan2(z_ecef, r_xy * (1.0 - e2 * N / (N + alt)));
    }

    double sin_lat = std::sin(lat);
    double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    double alt = r_xy / std::cos(lat) - N;

    return {
        lat * 180.0 / PI,
        lon * 180.0 / PI,
        alt
    };
}

int main() {
    std::cout << "=================================================\n";
    std::cout << "  SGP4 Propagator - ISS Tracking Example\n";
    std::cout << "=================================================\n\n";

    // ISS TLE - grab fresh ones from space-track.org or celestrak.org for better accuracy
    std::string tle_name = "ISS (ZARYA)";
    std::string tle_line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025";
    std::string tle_line2 = "2 25544  51.6400 208.9163 0006703 296.7361  63.2861 15.49560921426457";

    // Parse and initialize SGP4
    SGP4 sgp4;
    try {
        sgp4 = SGP4(tle_line1, tle_line2, tle_name);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing TLE: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Satellite: " << sgp4.tle.name << "\n";
    std::cout << "NORAD ID:  " << sgp4.tle.satellite_number << "\n";
    std::cout << "Epoch:     " << jdToDateString(sgp4.elems.epoch_jd) << "\n\n";

    // Display orbital parameters
    std::cout << "Orbital Elements:\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Inclination:     " << sgp4.tle.inclination << " deg\n";
    std::cout << "  RAAN:            " << sgp4.tle.raan << " deg\n";
    std::cout << "  Eccentricity:    " << std::setprecision(7) << sgp4.tle.eccentricity << "\n";
    std::cout << "  Arg of Perigee:  " << std::setprecision(4) << sgp4.tle.arg_perigee << " deg\n";
    std::cout << "  Mean Anomaly:    " << sgp4.tle.mean_anomaly << " deg\n";
    std::cout << "  Mean Motion:     " << std::setprecision(8) << sgp4.tle.mean_motion << " rev/day\n";
    std::cout << "  BSTAR:           " << std::scientific << sgp4.tle.bstar << std::fixed << "\n\n";

    std::cout << "Derived Parameters:\n";
    std::cout << std::setprecision(2);
    std::cout << "  Semi-major axis: " << sgp4.getSemiMajorAxis() << " km\n";
    std::cout << "  Orbital period:  " << sgp4.getPeriod() << " minutes\n";
    std::cout << "  Apogee altitude: " << sgp4.getApogee() << " km\n";
    std::cout << "  Perigee altitude:" << sgp4.getPerigee() << " km\n";
    std::cout << "  Deep space:      " << (sgp4.isDeepSpace() ? "Yes" : "No") << "\n\n";

    // =========================================================================
    // Propagate at epoch
    // =========================================================================
    std::cout << "State Vector at Epoch:\n";
    std::cout << "-----------------------\n";
    StateVector sv0 = sgp4.propagate(0.0);
    printStateVector(sv0);

    // =========================================================================
    // Propagate over one orbit
    // =========================================================================
    std::cout << "\nGround Track (1 orbit):\n";
    std::cout << "------------------------\n";
    std::cout << "Time [min]  Lat [deg]   Lon [deg]   Alt [km]\n";
    std::cout << "----------  ---------   ---------   --------\n";

    double period = sgp4.getPeriod();
    double dt = period / 20.0;  // 20 points per orbit

    for (int i = 0; i <= 20; ++i) {
        double t = i * dt;
        StateVector sv = sgp4.propagate(t);
        double jd = sgp4.elems.epoch_jd + t / constants::MIN_PER_DAY;
        GroundTrack gt = temeToGroundTrack(sv, jd);
        
        std::cout << std::setw(10) << std::setprecision(1) << t << "  "
                  << std::setw(9) << std::setprecision(2) << gt.lat << "   "
                  << std::setw(9) << gt.lon << "   "
                  << std::setw(8) << std::setprecision(1) << gt.alt << "\n";
    }

    // =========================================================================
    // Propagate over 24 hours
    // =========================================================================
    std::cout << "\n24-Hour Propagation Summary:\n";
    std::cout << "-----------------------------\n";

    std::vector<double> altitudes;
    std::vector<double> velocities;
    
    for (int i = 0; i <= 1440; ++i) {  // Every minute for 24 hours
        StateVector sv = sgp4.propagate(static_cast<double>(i));
        altitudes.push_back(sv.position_magnitude() - constants::RADIUS_EARTH);
        velocities.push_back(sv.velocity_magnitude());
    }

    auto minmax_alt = std::minmax_element(altitudes.begin(), altitudes.end());
    auto minmax_vel = std::minmax_element(velocities.begin(), velocities.end());

    std::cout << std::setprecision(2);
    std::cout << "  Altitude range: " << *minmax_alt.first << " - " << *minmax_alt.second << " km\n";
    std::cout << "  Velocity range: " << *minmax_vel.first << " - " << *minmax_vel.second << " km/s\n";
    std::cout << "  Orbits completed: " << std::setprecision(1) << 1440.0 / period << "\n";

    // =========================================================================
    // Compare states at different times
    // =========================================================================
    std::cout << "\nState Comparison at Different Times:\n";
    std::cout << "-------------------------------------\n";
    
    double times[] = {0.0, 30.0, 60.0, 90.0, 120.0};
    
    for (double t : times) {
        StateVector sv = sgp4.propagate(t);
        double jd = sgp4.elems.epoch_jd + t / constants::MIN_PER_DAY;
        GroundTrack gt = temeToGroundTrack(sv, jd);
        
        std::cout << "\nT + " << std::setprecision(0) << t << " min (" 
                  << jdToDateString(jd) << "):\n";
        std::cout << std::setprecision(3);
        std::cout << "  Position: [" << sv.x << ", " << sv.y << ", " << sv.z << "] km\n";
        std::cout << "  Velocity: [" << sv.vx << ", " << sv.vy << ", " << sv.vz << "] km/s\n";
        std::cout << std::setprecision(2);
        std::cout << "  Ground:   " << gt.lat << "° lat, " << gt.lon << "° lon, " 
                  << std::setprecision(1) << gt.alt << " km alt\n";
    }

    std::cout << "\n=================================================\n";
    std::cout << "  ISS Tracking complete!\n";
    std::cout << "=================================================\n";

    return 0;
}
