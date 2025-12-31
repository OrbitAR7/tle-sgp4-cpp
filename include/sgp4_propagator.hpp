#ifndef SGP4_PROPAGATOR_HPP
#define SGP4_PROPAGATOR_HPP

#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace sgp4 {

// ============================================================================
// Constants (WGS72 - used by SGP4)
// ============================================================================

namespace constants {
    // WGS72 Earth parameters (SGP4 uses WGS72, not WGS84)
    constexpr double MU = 398600.8;              // km^3/s^2
    constexpr double RADIUS_EARTH = 6378.135;    // km
    constexpr double J2 = 0.001082616;
    constexpr double J3 = -0.00000253881;
    constexpr double J4 = -0.00000165597;
    
    // Derived constants
    constexpr double XKE = 0.0743669161331734132;  // sqrt(MU) * (60/RADIUS_EARTH)^1.5
    constexpr double XKMPER = RADIUS_EARTH;
    constexpr double XMNPDA = 1440.0;              // Minutes per day
    constexpr double AE = 1.0;                     // Distance units/Earth radii
    constexpr double DE2RA = 0.0174532925199433;   // Degrees to radians
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWOPI = 2.0 * PI;
    constexpr double X3PIO2 = 3.0 * PI / 2.0;
    constexpr double PIO2 = PI / 2.0;
    
    // Time constants
    constexpr double MIN_PER_DAY = 1440.0;
    constexpr double SEC_PER_DAY = 86400.0;
    
    // Convergence
    constexpr double E6A = 1.0e-6;
}

// ============================================================================
// Data Structures
// ============================================================================

/**
 * @brief Position and velocity in TEME frame
 */
struct StateVector {
    double x, y, z;      // Position [km]
    double vx, vy, vz;   // Velocity [km/s]
    
    StateVector() : x(0), y(0), z(0), vx(0), vy(0), vz(0) {}
    StateVector(double x_, double y_, double z_, double vx_, double vy_, double vz_)
        : x(x_), y(y_), z(z_), vx(vx_), vy(vy_), vz(vz_) {}
    
    double position_magnitude() const { return std::sqrt(x*x + y*y + z*z); }
    double velocity_magnitude() const { return std::sqrt(vx*vx + vy*vy + vz*vz); }
};

/**
 * @brief Two-Line Element set data
 */
struct TLE {
    // Line 0 (optional name)
    std::string name;
    
    // Line 1 data
    int satellite_number;
    char classification;
    std::string intl_designator;
    int epoch_year;
    double epoch_day;
    double mean_motion_dot;      // First derivative of mean motion (rev/day^2)
    double mean_motion_ddot;     // Second derivative of mean motion (rev/day^3)
    double bstar;                // Drag term (1/Earth radii)
    int ephemeris_type;
    int element_set_number;
    
    // Line 2 data
    double inclination;          // [deg]
    double raan;                 // Right ascension of ascending node [deg]
    double eccentricity;         // [dimensionless]
    double arg_perigee;          // Argument of perigee [deg]
    double mean_anomaly;         // [deg]
    double mean_motion;          // [rev/day]
    int revolution_number;
    
    // Derived epoch in Julian Date
    double epoch_jd;
    
    TLE() : satellite_number(0), classification('U'), epoch_year(0), epoch_day(0),
            mean_motion_dot(0), mean_motion_ddot(0), bstar(0), ephemeris_type(0),
            element_set_number(0), inclination(0), raan(0), eccentricity(0),
            arg_perigee(0), mean_anomaly(0), mean_motion(0), revolution_number(0),
            epoch_jd(0) {}
};

/**
 * @brief SGP4 internal orbital elements (initialized from TLE)
 */
struct SGP4Elements {
    // Epoch
    double epoch_jd;
    
    // Mean elements at epoch (in radians and internal units)
    double no_kozai;    // Mean motion (rad/min) - Kozai mean motion
    double ecco;        // Eccentricity
    double inclo;       // Inclination (rad)
    double nodeo;       // RAAN (rad)
    double argpo;       // Argument of perigee (rad)
    double mo;          // Mean anomaly (rad)
    
    // Drag
    double bstar;
    
    // Derived quantities
    double no_unkozai;  // Mean motion without Kozai correction
    double a;           // Semi-major axis (Earth radii)
    double alta;        // Apogee altitude (Earth radii)
    double altp;        // Perigee altitude (Earth radii)
    
    // SGP4 initialization constants
    double d2, d3, d4;
    double t2cof, t3cof, t4cof, t5cof;
    double x1mth2, x7thm1;
    double xlcof, xmcof, xnodcf;
    double c1, c4, c5;
    double d2201, d2211, d3210, d3222, d4410, d4422;
    double d5220, d5232, d5421, d5433;
    double del1, del2, del3;
    double e3, ee2, se2, se3, sgh2, sgh3, sgh4;
    double sh2, sh3, si2, si3, sl2, sl3, sl4;
    double eta, etasq;
    double gsto;
    double omgcof;
    double xni, aycof;
    double xlmo, delmo, sinmao;
    
    // Deep space flag
    bool is_deep_space;
    bool is_simple;
    
    // Error flag
    int error;
    
    SGP4Elements() : epoch_jd(0), no_kozai(0), ecco(0), inclo(0), nodeo(0),
                     argpo(0), mo(0), bstar(0), no_unkozai(0), a(0), alta(0), altp(0),
                     d2(0), d3(0), d4(0), t2cof(0), t3cof(0), t4cof(0), t5cof(0),
                     x1mth2(0), x7thm1(0), xlcof(0), xmcof(0), xnodcf(0),
                     c1(0), c4(0), c5(0), d2201(0), d2211(0), d3210(0), d3222(0),
                     d4410(0), d4422(0), d5220(0), d5232(0), d5421(0), d5433(0),
                     del1(0), del2(0), del3(0), e3(0), ee2(0), se2(0), se3(0),
                     sgh2(0), sgh3(0), sgh4(0), sh2(0), sh3(0), si2(0), si3(0),
                     sl2(0), sl3(0), sl4(0), eta(0), etasq(0), gsto(0), omgcof(0),
                     xni(0), aycof(0), xlmo(0), delmo(0), sinmao(0),
                     is_deep_space(false), is_simple(false), error(0) {}
};

// ============================================================================
// TLE Parser
// ============================================================================

/**
 * @brief Parse exponential notation in TLE format (e.g., "12345-4" = 0.12345e-4)
 */
inline double parseTLEFloat(const std::string& s) {
    std::string str = s;
    // Remove leading/trailing whitespace
    size_t start = str.find_first_not_of(" ");
    size_t end = str.find_last_not_of(" ");
    if (start == std::string::npos) return 0.0;
    str = str.substr(start, end - start + 1);
    
    if (str.empty()) return 0.0;
    
    // Handle TLE exponential format: "12345-4" means 0.12345e-4
    // First check for standard format with 'e' or 'E'
    if (str.find('e') != std::string::npos || str.find('E') != std::string::npos) {
        return std::stod(str);
    }
    
    // Handle implicit decimal format
    size_t pos = str.find_first_of("+-", 1);  // Find sign after first char
    if (pos != std::string::npos && str[pos-1] != 'e' && str[pos-1] != 'E') {
        // Format: "12345-4" -> "0.12345e-4"
        std::string mantissa = str.substr(0, pos);
        std::string exponent = str.substr(pos);
        
        // Add leading "0." if not present
        if (mantissa[0] == '-' || mantissa[0] == '+') {
            return std::stod(mantissa.substr(0,1) + "0." + mantissa.substr(1) + "e" + exponent);
        } else {
            return std::stod("0." + mantissa + "e" + exponent);
        }
    }
    
    return std::stod(str);
}

/**
 * @brief Convert year and day-of-year to Julian Date
 */
inline double epochToJD(int year, double day) {
    // Handle two-digit year
    if (year < 57) {
        year += 2000;
    } else if (year < 100) {
        year += 1900;
    }
    
    // Julian date of Jan 1 of year
    int a = (14 - 1) / 12;
    int y = year + 4800 - a;
    int m = 1 + 12 * a - 3;
    
    double jd_jan1 = 1 + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    
    return jd_jan1 + day - 1.0;  // day is 1-based
}

/**
 * @brief Parse a Two-Line Element set from strings
 */
inline TLE parseTLE(const std::string& line1, const std::string& line2, 
                    const std::string& name = "") {
    TLE tle;
    tle.name = name;
    
    // Validate line lengths
    if (line1.length() < 69 || line2.length() < 69) {
        throw std::runtime_error("TLE lines too short");
    }
    
    // Validate line numbers
    if (line1[0] != '1' || line2[0] != '2') {
        throw std::runtime_error("Invalid TLE line numbers");
    }
    
    // Parse Line 1
    tle.satellite_number = std::stoi(line1.substr(2, 5));
    tle.classification = line1[7];
    tle.intl_designator = line1.substr(9, 8);
    
    // Epoch
    tle.epoch_year = std::stoi(line1.substr(18, 2));
    tle.epoch_day = std::stod(line1.substr(20, 12));
    
    // Mean motion derivatives
    tle.mean_motion_dot = std::stod(line1.substr(33, 10));
    tle.mean_motion_ddot = parseTLEFloat(line1.substr(44, 8));
    
    // BSTAR drag term
    tle.bstar = parseTLEFloat(line1.substr(53, 8));
    
    tle.ephemeris_type = std::stoi(line1.substr(62, 1));
    tle.element_set_number = std::stoi(line1.substr(64, 4));
    
    // Parse Line 2
    tle.inclination = std::stod(line2.substr(8, 8));
    tle.raan = std::stod(line2.substr(17, 8));
    
    // Eccentricity (implicit decimal point)
    std::string ecc_str = "0." + line2.substr(26, 7);
    tle.eccentricity = std::stod(ecc_str);
    
    tle.arg_perigee = std::stod(line2.substr(34, 8));
    tle.mean_anomaly = std::stod(line2.substr(43, 8));
    tle.mean_motion = std::stod(line2.substr(52, 11));
    tle.revolution_number = std::stoi(line2.substr(63, 5));
    
    // Compute Julian Date of epoch
    tle.epoch_jd = epochToJD(tle.epoch_year, tle.epoch_day);
    
    return tle;
}

// ============================================================================
// SGP4 Propagator
// ============================================================================

/**
 * @brief SGP4 orbit propagator
 */
class SGP4 {
public:
    SGP4Elements elems;
    TLE tle;
    
    SGP4() = default;
    
    /**
     * @brief Initialize SGP4 from TLE
     */
    explicit SGP4(const TLE& tle_) : tle(tle_) {
        initialize();
    }
    
    /**
     * @brief Initialize from TLE strings
     */
    SGP4(const std::string& line1, const std::string& line2, const std::string& name = "") {
        tle = parseTLE(line1, line2, name);
        initialize();
    }
    
    /**
     * @brief Propagate to minutes since epoch
     * @param tsince Minutes since TLE epoch
     * @return State vector in TEME frame
     */
    StateVector propagate(double tsince) const {
        return sgp4_propagate(tsince);
    }
    
    /**
     * @brief Propagate to Julian Date
     * @param jd Julian Date
     * @return State vector in TEME frame
     */
    StateVector propagateJD(double jd) const {
        double tsince = (jd - elems.epoch_jd) * constants::MIN_PER_DAY;
        return propagate(tsince);
    }
    
    /**
     * @brief Get semi-major axis [km]
     */
    double getSemiMajorAxis() const {
        return elems.a * constants::RADIUS_EARTH;
    }
    
    /**
     * @brief Get orbital period [minutes]
     */
    double getPeriod() const {
        return constants::TWOPI / elems.no_kozai;
    }
    
    /**
     * @brief Get apogee altitude [km]
     */
    double getApogee() const {
        return elems.alta * constants::RADIUS_EARTH;
    }
    
    /**
     * @brief Get perigee altitude [km]
     */
    double getPerigee() const {
        return elems.altp * constants::RADIUS_EARTH;
    }
    
    /**
     * @brief Check if deep space object
     */
    bool isDeepSpace() const {
        return elems.is_deep_space;
    }
    
    /**
     * @brief Greenwich Sidereal Time
     */
    static double gstime(double jd) {
        using namespace constants;
        
        double tut1 = (jd - 2451545.0) / 36525.0;
        double temp = -6.2e-6 * tut1 * tut1 * tut1 +
                      0.093104 * tut1 * tut1 +
                      (876600.0 * 3600.0 + 8640184.812866) * tut1 +
                      67310.54841;
        temp = std::fmod(temp * DE2RA / 240.0, TWOPI);
        
        if (temp < 0.0) temp += TWOPI;
        
        return temp;
    }
    
private:
    /**
     * @brief Initialize SGP4 elements from TLE
     */
    void initialize() {
        using namespace constants;
        
        elems.epoch_jd = tle.epoch_jd;
        elems.bstar = tle.bstar;
        
        // Convert to radians
        elems.inclo = tle.inclination * DE2RA;
        elems.nodeo = tle.raan * DE2RA;
        elems.argpo = tle.arg_perigee * DE2RA;
        elems.mo = tle.mean_anomaly * DE2RA;
        elems.ecco = tle.eccentricity;
        
        // Mean motion (rev/day -> rad/min)
        double no_raw = tle.mean_motion * TWOPI / MIN_PER_DAY;
        elems.no_kozai = no_raw;
        
        // Compute semi-major axis from mean motion
        double a1 = std::pow(XKE / no_raw, 2.0/3.0);
        double cosi = std::cos(elems.inclo);
        double d1 = 0.75 * J2 * (3.0 * cosi * cosi - 1.0) / 
                    (std::sqrt(1.0 - elems.ecco * elems.ecco) * 
                     (1.0 - elems.ecco * elems.ecco));
        double del_ = d1 / (a1 * a1);
        double adel = a1 * (1.0 - del_ * del_ - del_ * (1.0/3.0 + 134.0 * del_ * del_ / 81.0));
        del_ = d1 / (adel * adel);
        elems.no_unkozai = no_raw / (1.0 + del_);
        
        elems.a = std::pow(XKE / elems.no_unkozai, 2.0/3.0);
        elems.alta = elems.a * (1.0 + elems.ecco) - 1.0;
        elems.altp = elems.a * (1.0 - elems.ecco) - 1.0;
        
        // Check for deep space (period > 225 minutes)
        elems.is_deep_space = (TWOPI / elems.no_unkozai >= 225.0);
        
        if (elems.is_deep_space) {
            // Deep space initialization would go here
            // For simplicity, we'll handle near-Earth only in this implementation
            elems.error = 1;  // Mark as not fully supported
        }
        
        // Near-Earth initialization
        initNearEarth();
        
        // Compute Greenwich Sidereal Time at epoch
        elems.gsto = gstime(elems.epoch_jd);
    }
    
    /**
     * @brief Initialize near-Earth SGP4 constants
     */
    void initNearEarth() {
        using namespace constants;
        
        double cosio = std::cos(elems.inclo);
        double sinio = std::sin(elems.inclo);
        double cosio2 = cosio * cosio;
        
        double eccsq = elems.ecco * elems.ecco;
        double omeosq = 1.0 - eccsq;
        double rteosq = std::sqrt(omeosq);
        
        double posq = elems.a * omeosq;
        double rp = elems.a * (1.0 - elems.ecco);
        
        elems.x1mth2 = 1.0 - cosio2;
        elems.x7thm1 = 7.0 * cosio2 - 1.0;
        
        // Check for simple model (perigee < 220 km)
        elems.is_simple = (rp < (220.0 / XKMPER + 1.0));
        
        double con41 = -J2 * 1.5;
        double con42 = -J2 / 2.0;
        double x1m5th = 1.0 - 5.0 * cosio2;
        
        double ak = std::pow(XKE / elems.no_unkozai, 2.0/3.0);
        double d1_ = con41 / (ak * ak) / (posq * posq);
        
        double betao = std::sqrt(1.0 - eccsq);
        double betao2 = betao * betao;
        
        // Compute c1-c5
        double qoms24 = std::pow((120.0 - 78.0) / XKMPER, 4.0);  // (qo - s)^4
        double s4 = 78.0 / XKMPER;
        double perige = (rp - 1.0) * XKMPER;
        
        if (perige < 156.0) {
            s4 = perige - 78.0;
            if (perige < 98.0) s4 = 20.0;
            s4 = s4 / XKMPER + 1.0;
            qoms24 = std::pow((120.0 - s4 * XKMPER) / XKMPER, 4.0);
        }
        
        double pinvsq = 1.0 / (posq * posq);
        double tsi = 1.0 / (elems.a - s4);
        elems.eta = elems.a * elems.ecco * tsi;
        elems.etasq = elems.eta * elems.eta;
        double eeta = elems.ecco * elems.eta;
        double psisq = std::abs(1.0 - elems.etasq);
        double coef = qoms24 * std::pow(tsi, 4.0);
        double coef1 = coef / std::pow(psisq, 3.5);
        
        double c2 = coef1 * elems.no_unkozai * 
                    (elems.a * (1.0 + 1.5 * elems.etasq + eeta * (4.0 + elems.etasq)) +
                     0.375 * J2 * tsi / psisq * (3.0 * (1.0 - 3.0 * cosio2) * 
                     (1.0 + 1.5 * elems.etasq - 2.0 * eeta) +
                     0.75 * elems.x1mth2 * (2.0 * elems.etasq - eeta * (1.0 + elems.etasq)) *
                     std::cos(2.0 * elems.argpo)));
        
        elems.c1 = elems.bstar * c2;
        
        double c3 = 0.0;
        if (elems.ecco > 1.0e-4) {
            c3 = coef * tsi * J3 / J2 * elems.no_unkozai * sinio / elems.ecco;
        }
        
        elems.c4 = 2.0 * elems.no_unkozai * coef1 * elems.a * omeosq *
                   (elems.eta * (2.0 + 0.5 * elems.etasq) +
                    elems.ecco * (0.5 + 2.0 * elems.etasq) -
                    J2 * tsi / (elems.a * psisq) *
                    (-3.0 * (1.0 - 3.0 * cosio2) * (1.0 - 2.0 * eeta + elems.etasq * (1.5 - 0.5 * eeta)) +
                     0.75 * elems.x1mth2 * (2.0 * elems.etasq - eeta * (1.0 + elems.etasq)) *
                     std::cos(2.0 * elems.argpo)));
        
        elems.c5 = 2.0 * coef1 * elems.a * omeosq *
                   (1.0 + 2.75 * (elems.etasq + eeta) + eeta * elems.etasq);
        
        // Additional coefficients
        double theta4 = cosio2 * cosio2;
        double temp1 = 3.0 * J2 * pinvsq * elems.no_unkozai;
        double temp2 = temp1 * J2 * pinvsq;
        double temp3 = 1.25 * J4 * pinvsq * pinvsq * elems.no_unkozai;
        
        elems.xmcof = 0.0;
        if (elems.ecco > 1.0e-4) {
            elems.xmcof = -2.0/3.0 * coef * elems.bstar / eeta;
        }
        
        elems.xlcof = 0.125 * J3 / J2 * sinio * 
                      (3.0 + 5.0 * cosio) / (1.0 + cosio);
        elems.aycof = 0.25 * J3 / J2 * sinio;
        elems.xnodcf = 3.5 * omeosq * J2 * pinvsq * elems.no_unkozai;
        
        elems.t2cof = 1.5 * elems.c1;
        
        if (!elems.is_simple) {
            double d2_ = 4.0 * elems.a * tsi * elems.c1 * elems.c1;
            double temp = d2_ * tsi * elems.c1 / 3.0;
            elems.d2 = d2_;
            elems.d3 = (17.0 * elems.a + s4) * temp;
            elems.d4 = 0.5 * temp * elems.a * tsi * 
                       (221.0 * elems.a + 31.0 * s4) * elems.c1;
            elems.t3cof = elems.d2 + 2.0 * elems.c1 * elems.c1;
            elems.t4cof = 0.25 * (3.0 * elems.d3 + elems.c1 * 
                          (12.0 * elems.d2 + 10.0 * elems.c1 * elems.c1));
            elems.t5cof = 0.2 * (3.0 * elems.d4 + 12.0 * elems.c1 * elems.d3 +
                          6.0 * elems.d2 * elems.d2 + 
                          15.0 * elems.c1 * elems.c1 * (2.0 * elems.d2 + elems.c1 * elems.c1));
        }
        
        // Store for mean anomaly update
        elems.delmo = std::pow(1.0 + elems.eta * std::cos(elems.mo), 3.0);
        elems.sinmao = std::sin(elems.mo);
    }
    
    /**
     * @brief SGP4 propagation
     */
    StateVector sgp4_propagate(double tsince) const {
        using namespace constants;
        
        StateVector sv;
        
        // Update for secular effects
        double xmdf = elems.mo + elems.no_kozai * tsince;
        double argpdf = elems.argpo;
        double nodedf = elems.nodeo;
        
        double cosio = std::cos(elems.inclo);
        double sinio = std::sin(elems.inclo);
        
        double t2 = tsince * tsince;
        double tempa = 1.0 - elems.c1 * tsince;
        double tempe = elems.bstar * elems.c4 * tsince;
        double templ = elems.t2cof * t2;
        
        if (!elems.is_simple) {
            double t3 = t2 * tsince;
            double t4 = t3 * tsince;
            tempa -= elems.d2 * t2 - elems.d3 * t3 - elems.d4 * t4;
            tempe += elems.bstar * elems.c5 * (std::sin(xmdf) - elems.sinmao);
            templ += elems.t3cof * t3 + t4 * (elems.t4cof + tsince * elems.t5cof);
        }
        
        double a = elems.a * tempa * tempa;
        double e = elems.ecco - tempe;
        double xl = xmdf + argpdf + nodedf + elems.no_kozai * templ;
        
        // Fix tolerance
        if (e < 1.0e-6) e = 1.0e-6;
        if (e > 0.999) e = 0.999;
        
        double beta = std::sqrt(1.0 - e * e);
        double n = XKE / std::pow(a, 1.5);
        
        // Long period periodics
        double axn = e * std::cos(argpdf);
        double ayn = e * std::sin(argpdf) + elems.aycof;
        double xlt = xl + elems.xlcof * axn;
        
        // Solve Kepler's equation
        double u = std::fmod(xlt - nodedf, TWOPI);
        double eo1 = u;
        double tem5 = 1.0;
        int iter = 0;
        
        while (std::abs(tem5) >= E6A && iter < 10) {
            double sineo1 = std::sin(eo1);
            double coseo1 = std::cos(eo1);
            tem5 = 1.0 - coseo1 * axn - sineo1 * ayn;
            tem5 = (u - ayn * coseo1 + axn * sineo1 - eo1) / tem5;
            eo1 += tem5;
            iter++;
        }
        
        double sineo1 = std::sin(eo1);
        double coseo1 = std::cos(eo1);
        
        // Short period preliminary quantities
        double ecose = axn * coseo1 + ayn * sineo1;
        double esine = axn * sineo1 - ayn * coseo1;
        double el2 = axn * axn + ayn * ayn;
        double pl = a * (1.0 - el2);
        
        if (pl < 0.0) {
            sv = StateVector();
            return sv;
        }
        
        double r = a * (1.0 - ecose);
        double rdot = XKE * std::sqrt(a) * esine / r;
        double rvdot = XKE * std::sqrt(pl) / r;
        double betal = std::sqrt(1.0 - el2);
        
        double sinu = a / r * (sineo1 - ayn - axn * esine / (1.0 + betal));
        double cosu = a / r * (coseo1 - axn + ayn * esine / (1.0 + betal));
        double su = std::atan2(sinu, cosu);
        
        double sin2u = (cosu + cosu) * sinu;
        double cos2u = 1.0 - 2.0 * sinu * sinu;
        
        // Short periodics
        double rk = r * (1.0 - 1.5 * J2 * std::sqrt(1.0 - elems.ecco * elems.ecco) / 
                   (pl * pl) * (3.0 * cosio * cosio - 1.0)) +
                   0.5 * J2 / (pl * pl) * elems.x1mth2 * cos2u;
        double uk = su - 0.25 * J2 / (pl * pl) * elems.x7thm1 * sin2u;
        double xnodek = nodedf + 1.5 * J2 / (pl * pl) * cosio * sin2u;
        double xinck = elems.inclo + 1.5 * J2 / (pl * pl) * cosio * sinio * cos2u;
        double rdotk = rdot - J2 * n / (pl) * elems.x1mth2 * sin2u;
        double rfdotk = rvdot + J2 * n / (pl) * (elems.x1mth2 * cos2u + 
                        1.5 * (1.0 - 3.0 * cosio * cosio));
        
        // Orientation vectors
        double sinuk = std::sin(uk);
        double cosuk = std::cos(uk);
        double sinik = std::sin(xinck);
        double cosik = std::cos(xinck);
        double sinnok = std::sin(xnodek);
        double cosnok = std::cos(xnodek);
        
        double xmx = -sinnok * cosik;
        double xmy = cosnok * cosik;
        double ux = xmx * sinuk + cosnok * cosuk;
        double uy = xmy * sinuk + sinnok * cosuk;
        double uz = sinik * sinuk;
        double vx = xmx * cosuk - cosnok * sinuk;
        double vy = xmy * cosuk - sinnok * sinuk;
        double vz = sinik * cosuk;
        
        // Position and velocity in km and km/s
        sv.x = rk * ux * XKMPER;
        sv.y = rk * uy * XKMPER;
        sv.z = rk * uz * XKMPER;
        sv.vx = (rdotk * ux + rfdotk * vx) * XKMPER / 60.0;
        sv.vy = (rdotk * uy + rfdotk * vy) * XKMPER / 60.0;
        sv.vz = (rdotk * uz + rfdotk * vz) * XKMPER / 60.0;
        
        return sv;
    }
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert Julian Date to calendar date string
 */
inline std::string jdToDateString(double jd) {
    // Algorithm from Meeus
    double z = std::floor(jd + 0.5);
    double f = jd + 0.5 - z;
    
    double a;
    if (z < 2299161) {
        a = z;
    } else {
        double alpha = std::floor((z - 1867216.25) / 36524.25);
        a = z + 1 + alpha - std::floor(alpha / 4.0);
    }
    
    double b = a + 1524;
    double c = std::floor((b - 122.1) / 365.25);
    double d = std::floor(365.25 * c);
    double e = std::floor((b - d) / 30.6001);
    
    int day = static_cast<int>(b - d - std::floor(30.6001 * e));
    int month = static_cast<int>((e < 14) ? e - 1 : e - 13);
    int year = static_cast<int>((month > 2) ? c - 4716 : c - 4715);
    
    double hours_dec = f * 24.0;
    int hours = static_cast<int>(hours_dec);
    int minutes = static_cast<int>((hours_dec - hours) * 60);
    int seconds = static_cast<int>(((hours_dec - hours) * 60 - minutes) * 60);
    
    std::ostringstream oss;
    oss << year << "-" << std::setfill('0') << std::setw(2) << month << "-"
        << std::setw(2) << day << " " << std::setw(2) << hours << ":"
        << std::setw(2) << minutes << ":" << std::setw(2) << seconds << " UTC";
    
    return oss.str();
}

/**
 * @brief Print state vector
 */
inline void printStateVector(const StateVector& sv, const std::string& label = "") {
    if (!label.empty()) {
        std::cout << label << std::endl;
    }
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Position [km]: " << sv.x << ", " << sv.y << ", " << sv.z << std::endl;
    std::cout << "  Velocity [km/s]: " << sv.vx << ", " << sv.vy << ", " << sv.vz << std::endl;
    std::cout << "  |r| = " << sv.position_magnitude() << " km" << std::endl;
    std::cout << "  |v| = " << sv.velocity_magnitude() << " km/s" << std::endl;
}

} // namespace sgp4

#endif // SGP4_PROPAGATOR_HPP
