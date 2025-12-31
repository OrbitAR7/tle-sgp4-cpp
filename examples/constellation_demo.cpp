/**
 * Multi-Satellite Tracking Demo
 * 
 * Shows how to track multiple satellites at once - LEO, MEO, and HEO orbits.
 * Pretty cool to see how different orbit types behave over 24 hours.
 * 
 * Author: Hamoud Alshammari
 * 2025
 */

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include "../include/sgp4_propagator.hpp"

using namespace sgp4;

// Sample TLEs for various satellites
struct TLEData {
    std::string name;
    std::string line1;
    std::string line2;
};

int main() {
    std::cout << "=================================================\n";
    std::cout << "  SGP4 Propagator - Multi-Satellite Demo\n";
    std::cout << "=================================================\n\n";

    // Sample TLE set - mix of different orbit types to show the variety
    std::vector<TLEData> tle_database = {
        {
            "ISS (ZARYA)",
            "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025",
            "2 25544  51.6400 208.9163 0006703 296.7361  63.2861 15.49560921426457"
        },
        {
            "HUBBLE",
            "1 20580U 90037B   24001.50000000  .00000800  00000-0  40000-4 0  9999",
            "2 20580  28.4700 120.5000 0002800 100.0000 260.0000 15.09000000100000"
        },
        {
            "GPS BIIR-2",
            "1 24876U 97035A   24001.50000000 -.00000020  00000-0  00000+0 0  9999",
            "2 24876  55.7000  50.0000 0100000 270.0000  90.0000  2.00560000100000"
        },
        {
            "MOLNIYA 1-91",
            "1 25485U 98054A   24001.50000000  .00000100  00000-0  10000-4 0  9999",
            "2 25485  62.8000 280.0000 7000000 270.0000  90.0000  2.00600000100000"
        },
        {
            "STARLINK-1234",
            "1 45000U 20001A   24001.50000000  .00001000  00000-0  50000-4 0  9999",
            "2 45000  53.0000 100.0000 0001500  90.0000 270.0000 15.40000000100000"
        }
    };

    // Parse all TLEs
    std::vector<SGP4> satellites;
    
    std::cout << "Loading satellites...\n\n";
    
    for (const auto& tle_data : tle_database) {
        try {
            SGP4 sat(tle_data.line1, tle_data.line2, tle_data.name);
            satellites.push_back(sat);
            std::cout << "  ✓ " << tle_data.name << " (NORAD " << sat.tle.satellite_number << ")\n";
        } catch (const std::exception& e) {
            std::cout << "  ✗ " << tle_data.name << " - Error: " << e.what() << "\n";
        }
    }

    std::cout << "\nLoaded " << satellites.size() << " satellites.\n\n";

    // =========================================================================
    // Compare orbital characteristics
    // =========================================================================
    std::cout << "Orbital Characteristics:\n";
    std::cout << std::string(100, '-') << "\n";
    std::cout << std::left << std::setw(20) << "Satellite"
              << std::right << std::setw(12) << "Inc [deg]"
              << std::setw(12) << "Ecc"
              << std::setw(14) << "Period [min]"
              << std::setw(14) << "Apogee [km]"
              << std::setw(14) << "Perigee [km]"
              << std::setw(10) << "Type"
              << "\n";
    std::cout << std::string(100, '-') << "\n";

    for (const auto& sat : satellites) {
        std::string orbit_type;
        double period = sat.getPeriod();
        
        // Classify with Molniya/HEO preference before period buckets
        if (sat.tle.eccentricity > 0.5 && std::abs(sat.tle.inclination - 63.4) < 5.0) {
            orbit_type = "Molniya";  // high-ecc, critical inclination
        } else if (sat.tle.eccentricity > 0.1) {
            orbit_type = "HEO";
        } else if (period > 1400) {
            orbit_type = "GEO";
        } else if (period > 600) {
            orbit_type = "MEO";
        } else {
            orbit_type = "LEO";
        }
        
        std::cout << std::left << std::setw(20) << sat.tle.name
                  << std::right << std::fixed
                  << std::setw(12) << std::setprecision(2) << sat.tle.inclination
                  << std::setw(12) << std::setprecision(6) << sat.tle.eccentricity
                  << std::setw(14) << std::setprecision(2) << period
                  << std::setw(14) << std::setprecision(1) << sat.getApogee()
                  << std::setw(14) << sat.getPerigee()
                  << std::setw(10) << orbit_type
                  << "\n";
    }

    std::cout << std::string(100, '-') << "\n\n";

    // =========================================================================
    // Propagate all satellites to a common time
    // =========================================================================
    std::cout << "State Vectors at Epoch:\n";
    std::cout << std::string(80, '-') << "\n";

    for (const auto& sat : satellites) {
        StateVector sv = sat.propagate(0.0);
        
        std::cout << sat.tle.name << ":\n";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  r = [" << std::setw(12) << sv.x << ", " 
                  << std::setw(12) << sv.y << ", " 
                  << std::setw(12) << sv.z << "] km\n";
        std::cout << "  v = [" << std::setw(12) << sv.vx << ", " 
                  << std::setw(12) << sv.vy << ", " 
                  << std::setw(12) << sv.vz << "] km/s\n";
        std::cout << "  |r| = " << std::setprecision(1) << sv.position_magnitude() 
                  << " km, |v| = " << std::setprecision(3) << sv.velocity_magnitude() << " km/s\n\n";
    }

    // =========================================================================
    // Time evolution over 24 hours
    // =========================================================================
    std::cout << "Altitude Evolution Over 24 Hours:\n";
    std::cout << std::string(80, '-') << "\n";

    // Machine-readable header for plotting
    std::cout << "#ALT_NAMES";
    for (const auto& sat : satellites) {
        std::string n = sat.tle.name;
        std::replace(n.begin(), n.end(), ' ', '_');
        std::cout << "," << n;
    }
    std::cout << "\n";
    
    std::cout << std::left << std::setw(8) << "T [min]";
    for (const auto& sat : satellites) {
        // Truncate name if too long
        std::string short_name = sat.tle.name.substr(0, 12);
        std::cout << std::right << std::setw(14) << short_name;
    }
    std::cout << "\n";
    std::cout << std::string(8 + satellites.size() * 14, '-') << "\n";

    const int horizon_minutes = 24 * 60;  // 1 day
    const int step_minutes = 30;           // sampling step
    for (int t = 0; t <= horizon_minutes; t += step_minutes) {
        std::cout << std::left << std::setw(8) << t;
        
        for (const auto& sat : satellites) {
            StateVector sv = sat.propagate(static_cast<double>(t));
            double alt = sv.position_magnitude() - constants::RADIUS_EARTH;
            std::cout << std::right << std::fixed << std::setprecision(1) 
                      << std::setw(14) << alt;
        }
        std::cout << "\n";
    }

    // =========================================================================
    // Find closest approach between ISS and other satellites
    // =========================================================================
    if (satellites.size() >= 2) {
        std::cout << "\nClosest Approach Analysis (ISS vs others, 24 hours):\n";
        std::cout << std::string(60, '-') << "\n";
        
        const SGP4& iss = satellites[0];  // Assuming ISS is first
        
        for (size_t i = 1; i < satellites.size(); ++i) {
            const SGP4& other = satellites[i];
            
            double min_dist = 1e9;
            double min_time = 0;
            
            // Check every minute for 24 hours
            for (int t = 0; t <= 1440; ++t) {
                StateVector sv_iss = iss.propagate(static_cast<double>(t));
                StateVector sv_other = other.propagate(static_cast<double>(t));
                
                double dx = sv_iss.x - sv_other.x;
                double dy = sv_iss.y - sv_other.y;
                double dz = sv_iss.z - sv_other.z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist < min_dist) {
                    min_dist = dist;
                    min_time = t;
                }
            }
            
            std::cout << std::left << std::setw(20) << other.tle.name
                      << " Min distance: " << std::right << std::fixed 
                      << std::setprecision(1) << std::setw(10) << min_dist << " km"
                      << " at T+" << std::setw(6) << min_time << " min\n";
        }
    }

    std::cout << "\n=================================================\n";
    std::cout << "  Multi-Satellite Demo complete!\n";
    std::cout << "=================================================\n";

    return 0;
}
