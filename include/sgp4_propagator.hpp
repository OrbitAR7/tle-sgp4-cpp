#ifndef SGP4_PROPAGATOR_HPP
#define SGP4_PROPAGATOR_HPP

// Near-Earth SGP4 (orbital period < 225 min), following
// Vallado, Crawford, Hujsak & Kelso, "Revisiting Spacetrack Report #3",
// AIAA 2006-6753, with WGS-72 constants and the "improved" operation mode.
// Deep-space objects (SDP4) are detected and rejected, not propagated.

#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace sgp4 {

// ============================================================================
// Constants (WGS-72, as used by SGP4)
// ============================================================================

namespace constants {
    constexpr double MU = 398600.8;              // km^3/s^2
    constexpr double RADIUS_EARTH = 6378.135;    // km
    constexpr double J2 = 0.001082616;
    constexpr double J3 = -0.00000253881;
    constexpr double J4 = -0.00000165597;
    constexpr double J3OJ2 = J3 / J2;

    // sqrt(MU / RADIUS_EARTH^3) in 1/min, i.e. 60 / sqrt(RE^3 / MU)
    constexpr double XKE = 0.0743669161331734132;
    constexpr double XKMPER = RADIUS_EARTH;
    constexpr double VKMPERSEC = RADIUS_EARTH * XKE / 60.0;
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWOPI = 2.0 * PI;
    constexpr double DE2RA = PI / 180.0;
    constexpr double X2O3 = 2.0 / 3.0;

    constexpr double MIN_PER_DAY = 1440.0;
    constexpr double SEC_PER_DAY = 86400.0;
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
 * @brief SGP4 elements and initialization constants (from sgp4init)
 */
struct SGP4Elements {
    double epoch_jd = 0;

    // Mean elements at epoch (radians, rad/min)
    double no_kozai = 0, no_unkozai = 0;
    double ecco = 0, inclo = 0, nodeo = 0, argpo = 0, mo = 0;
    double bstar = 0;

    double a = 0, alta = 0, altp = 0;      // Earth radii
    double con41 = 0, x1mth2 = 0, x7thm1 = 0;
    double cc1 = 0, cc4 = 0, cc5 = 0;
    double d2 = 0, d3 = 0, d4 = 0;
    double t2cof = 0, t3cof = 0, t4cof = 0, t5cof = 0;
    double mdot = 0, argpdot = 0, nodedot = 0, nodecf = 0;
    double omgcof = 0, xmcof = 0, xlcof = 0, aycof = 0;
    double eta = 0, delmo = 0, sinmao = 0;
    double gsto = 0;

    bool is_deep_space = false;
    bool is_simple = false;   // perigee < 220 km: drop the higher-order drag terms
    int error = 0;            // 0 ok, 1 deep space (unsupported)
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
    
    // Julian day number of Jan 1 (an integer, referring to noon)
    int jdn_jan1 = 1 + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    
    // JD at 0h on Jan 1 is JDN - 0.5; TLE day-of-year is 1-based
    return jdn_jan1 - 0.5 + day - 1.0;
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
 * @brief Near-Earth SGP4 orbit propagator. Output is in the TEME frame.
 */
class SGP4 {
public:
    SGP4Elements elems;
    TLE tle;

    SGP4() = default;

    explicit SGP4(const TLE& tle_) : tle(tle_) {
        initialize();
    }

    SGP4(const std::string& line1, const std::string& line2, const std::string& name = "") {
        tle = parseTLE(line1, line2, name);
        initialize();
    }

    /**
     * @brief Propagate to minutes since epoch
     * @throws std::domain_error for deep-space objects or a decayed/invalid orbit
     */
    StateVector propagate(double tsince) const {
        if (elems.is_deep_space) {
            throw std::domain_error("SGP4: deep-space object (period >= 225 min); SDP4 is not implemented");
        }
        return sgp4_propagate(tsince);
    }

    StateVector propagateJD(double jd) const {
        return propagate((jd - elems.epoch_jd) * constants::MIN_PER_DAY);
    }

    double getSemiMajorAxis() const { return elems.a * constants::RADIUS_EARTH; }        // km
    double getPeriod() const { return constants::TWOPI / elems.no_unkozai; }             // min
    double getApogee() const { return elems.alta * constants::RADIUS_EARTH; }            // km altitude
    double getPerigee() const { return elems.altp * constants::RADIUS_EARTH; }           // km altitude
    bool isDeepSpace() const { return elems.is_deep_space; }

    /**
     * @brief Greenwich mean sidereal time [rad] (IAU-82, UT1 ~ UTC)
     */
    static double gstime(double jdut1) {
        using namespace constants;
        double tut1 = (jdut1 - 2451545.0) / 36525.0;
        double temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
                      (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;  // sec
        temp = std::fmod(temp * DE2RA / 240.0, TWOPI);
        if (temp < 0.0) temp += TWOPI;
        return temp;
    }

private:
    void initialize() {
        using namespace constants;
        SGP4Elements& e = elems;

        e.epoch_jd = tle.epoch_jd;
        e.bstar = tle.bstar;
        e.inclo = tle.inclination * DE2RA;
        e.nodeo = tle.raan * DE2RA;
        e.argpo = tle.arg_perigee * DE2RA;
        e.mo = tle.mean_anomaly * DE2RA;
        e.ecco = tle.eccentricity;
        e.no_kozai = tle.mean_motion * TWOPI / MIN_PER_DAY;   // rad/min

        // ---- initl: recover Brouwer mean motion and semi-major axis
        const double eccsq = e.ecco * e.ecco;
        const double omeosq = 1.0 - eccsq;
        const double rteosq = std::sqrt(omeosq);
        const double cosio = std::cos(e.inclo);
        const double cosio2 = cosio * cosio;
        const double sinio = std::sin(e.inclo);

        const double ak = std::pow(XKE / e.no_kozai, X2O3);
        const double d1 = 0.75 * J2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
        double del = d1 / (ak * ak);
        const double adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
        del = d1 / (adel * adel);
        e.no_unkozai = e.no_kozai / (1.0 + del);

        const double ao = std::pow(XKE / e.no_unkozai, X2O3);
        const double po = ao * omeosq;
        const double con42 = 1.0 - 5.0 * cosio2;
        e.con41 = -con42 - cosio2 - cosio2;          // 3 cos^2 i - 1
        const double posq = po * po;
        const double rp = ao * (1.0 - e.ecco);

        e.a = ao;
        e.alta = ao * (1.0 + e.ecco) - 1.0;
        e.altp = rp - 1.0;
        e.gsto = gstime(e.epoch_jd);

        e.is_deep_space = (TWOPI / e.no_unkozai >= 225.0);
        if (e.is_deep_space) {
            e.error = 1;
            return;
        }

        // ---- sgp4init (near-Earth)
        e.is_simple = rp < (220.0 / RADIUS_EARTH + 1.0);

        double sfour = 78.0 / RADIUS_EARTH + 1.0;
        double qzms24 = std::pow((120.0 - 78.0) / RADIUS_EARTH, 4.0);
        const double perige = (rp - 1.0) * RADIUS_EARTH;
        if (perige < 156.0) {
            sfour = perige - 78.0;
            if (perige < 98.0) sfour = 20.0;
            qzms24 = std::pow((120.0 - sfour) / RADIUS_EARTH, 4.0);
            sfour = sfour / RADIUS_EARTH + 1.0;
        }

        const double pinvsq = 1.0 / posq;
        const double tsi = 1.0 / (ao - sfour);
        e.eta = ao * e.ecco * tsi;
        const double etasq = e.eta * e.eta;
        const double eeta = e.ecco * e.eta;
        const double psisq = std::fabs(1.0 - etasq);
        const double coef = qzms24 * std::pow(tsi, 4.0);
        const double coef1 = coef / std::pow(psisq, 3.5);

        const double cc2 = coef1 * e.no_unkozai *
            (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
             0.375 * J2 * tsi / psisq * e.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
        e.cc1 = e.bstar * cc2;
        double cc3 = 0.0;
        if (e.ecco > 1.0e-4) {
            cc3 = -2.0 * coef * tsi * J3OJ2 * e.no_unkozai * sinio / e.ecco;
        }
        e.x1mth2 = 1.0 - cosio2;
        e.cc4 = 2.0 * e.no_unkozai * coef1 * ao * omeosq *
            (e.eta * (2.0 + 0.5 * etasq) + e.ecco * (0.5 + 2.0 * etasq) -
             J2 * tsi / (ao * psisq) *
                 (-3.0 * e.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                  0.75 * e.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * std::cos(2.0 * e.argpo)));
        e.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

        // Secular rates from J2 and J4
        const double cosio4 = cosio2 * cosio2;
        const double temp1 = 1.5 * J2 * pinvsq * e.no_unkozai;
        const double temp2 = 0.5 * temp1 * J2 * pinvsq;
        const double temp3 = -0.46875 * J4 * pinvsq * pinvsq * e.no_unkozai;
        e.mdot = e.no_unkozai + 0.5 * temp1 * rteosq * e.con41 +
                 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
        e.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                    temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
        const double xhdot1 = -temp1 * cosio;
        e.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;

        e.omgcof = e.bstar * cc3 * std::cos(e.argpo);
        e.xmcof = (e.ecco > 1.0e-4) ? -X2O3 * coef * e.bstar / eeta : 0.0;
        e.nodecf = 3.5 * omeosq * xhdot1 * e.cc1;
        e.t2cof = 1.5 * e.cc1;
        const double denom = (std::fabs(cosio + 1.0) > 1.5e-12) ? (1.0 + cosio) : 1.5e-12;
        e.xlcof = -0.25 * J3OJ2 * sinio * (3.0 + 5.0 * cosio) / denom;
        e.aycof = -0.5 * J3OJ2 * sinio;
        e.delmo = std::pow(1.0 + e.eta * std::cos(e.mo), 3.0);
        e.sinmao = std::sin(e.mo);
        e.x7thm1 = 7.0 * cosio2 - 1.0;

        if (!e.is_simple) {
            const double cc1sq = e.cc1 * e.cc1;
            e.d2 = 4.0 * ao * tsi * cc1sq;
            const double temp = e.d2 * tsi * e.cc1 / 3.0;
            e.d3 = (17.0 * ao + sfour) * temp;
            e.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * e.cc1;
            e.t3cof = e.d2 + 2.0 * cc1sq;
            e.t4cof = 0.25 * (3.0 * e.d3 + e.cc1 * (12.0 * e.d2 + 10.0 * cc1sq));
            e.t5cof = 0.2 * (3.0 * e.d4 + 12.0 * e.cc1 * e.d3 + 6.0 * e.d2 * e.d2 +
                             15.0 * cc1sq * (2.0 * e.d2 + cc1sq));
        }
    }

    StateVector sgp4_propagate(double t) const {
        using namespace constants;
        const SGP4Elements& e = elems;

        // ---- secular gravity and atmospheric drag
        const double xmdf = e.mo + e.mdot * t;
        const double argpdf = e.argpo + e.argpdot * t;
        const double nodedf = e.nodeo + e.nodedot * t;
        double argpm = argpdf;
        double mm = xmdf;
        const double t2 = t * t;
        double nodem = nodedf + e.nodecf * t2;
        double tempa = 1.0 - e.cc1 * t;
        double tempe = e.bstar * e.cc4 * t;
        double templ = e.t2cof * t2;

        if (!e.is_simple) {
            const double delomg = e.omgcof * t;
            const double delmtemp = 1.0 + e.eta * std::cos(xmdf);
            const double delm = e.xmcof * (delmtemp * delmtemp * delmtemp - e.delmo);
            const double temp = delomg + delm;
            mm = xmdf + temp;
            argpm = argpdf - temp;
            const double t3 = t2 * t;
            const double t4 = t3 * t;
            tempa = tempa - e.d2 * t2 - e.d3 * t3 - e.d4 * t4;
            tempe = tempe + e.bstar * e.cc5 * (std::sin(mm) - e.sinmao);
            templ = templ + e.t3cof * t3 + t4 * (e.t4cof + t * e.t5cof);
        }

        double nm = e.no_unkozai;
        double em = e.ecco;
        const double inclm = e.inclo;
        const double am = std::pow(XKE / nm, X2O3) * tempa * tempa;
        nm = XKE / std::pow(am, 1.5);
        em = em - tempe;
        if (em >= 1.0 || em < -0.001 || am < 0.95) {
            throw std::domain_error("SGP4: eccentricity/semi-major axis out of range (decayed orbit?)");
        }
        if (em < 1.0e-6) em = 1.0e-6;
        mm = mm + e.no_unkozai * templ;
        double xlm = mm + argpm + nodem;

        nodem = std::fmod(nodem, TWOPI);
        argpm = std::fmod(argpm, TWOPI);
        xlm = std::fmod(xlm, TWOPI);
        mm = std::fmod(xlm - argpm - nodem, TWOPI);

        const double sinip = std::sin(inclm);
        const double cosip = std::cos(inclm);

        // ---- long-period periodics
        const double axnl = em * std::cos(argpm);
        double temp = 1.0 / (am * (1.0 - em * em));
        const double aynl = em * std::sin(argpm) + temp * e.aycof;
        const double xl = mm + argpm + nodem + temp * e.xlcof * axnl;

        // ---- Kepler's equation
        const double u = std::fmod(xl - nodem, TWOPI);
        double eo1 = u;
        double tem5 = 9999.9;
        double sineo1 = 0.0, coseo1 = 0.0;
        for (int ktr = 1; std::fabs(tem5) >= 1.0e-12 && ktr <= 10; ++ktr) {
            sineo1 = std::sin(eo1);
            coseo1 = std::cos(eo1);
            tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
            if (std::fabs(tem5) >= 0.95) tem5 = tem5 > 0.0 ? 0.95 : -0.95;
            eo1 = eo1 + tem5;
        }

        // ---- short-period preliminaries
        const double ecose = axnl * coseo1 + aynl * sineo1;
        const double esine = axnl * sineo1 - aynl * coseo1;
        const double el2 = axnl * axnl + aynl * aynl;
        const double pl = am * (1.0 - el2);
        if (pl < 0.0) {
            throw std::domain_error("SGP4: semi-latus rectum < 0");
        }
        const double rl = am * (1.0 - ecose);
        const double rdotl = std::sqrt(am) * esine / rl;
        const double rvdotl = std::sqrt(pl) / rl;
        const double betal = std::sqrt(1.0 - el2);
        temp = esine / (1.0 + betal);
        const double sinu = am / rl * (sineo1 - aynl - axnl * temp);
        const double cosu = am / rl * (coseo1 - axnl + aynl * temp);
        double su = std::atan2(sinu, cosu);
        const double sin2u = (cosu + cosu) * sinu;
        const double cos2u = 1.0 - 2.0 * sinu * sinu;
        temp = 1.0 / pl;
        const double temp1 = 0.5 * J2 * temp;
        const double temp2 = temp1 * temp;

        // ---- short-period periodics
        const double mrt = rl * (1.0 - 1.5 * temp2 * betal * e.con41) + 0.5 * temp1 * e.x1mth2 * cos2u;
        su = su - 0.25 * temp2 * e.x7thm1 * sin2u;
        const double xnode = nodem + 1.5 * temp2 * cosip * sin2u;
        const double xinc = inclm + 1.5 * temp2 * cosip * sinip * cos2u;
        const double mvt = rdotl - nm * temp1 * e.x1mth2 * sin2u / XKE;
        const double rvdot = rvdotl + nm * temp1 * (e.x1mth2 * cos2u + 1.5 * e.con41) / XKE;

        // ---- orientation vectors
        const double sinsu = std::sin(su), cossu = std::cos(su);
        const double snod = std::sin(xnode), cnod = std::cos(xnode);
        const double sini = std::sin(xinc), cosi = std::cos(xinc);
        const double xmx = -snod * cosi;
        const double xmy = cnod * cosi;
        const double ux = xmx * sinsu + cnod * cossu;
        const double uy = xmy * sinsu + snod * cossu;
        const double uz = sini * sinsu;
        const double vx = xmx * cossu - cnod * sinsu;
        const double vy = xmy * cossu - snod * sinsu;
        const double vz = sini * cossu;

        if (mrt < 1.0) {
            throw std::domain_error("SGP4: satellite has decayed");
        }

        return StateVector(mrt * ux * XKMPER, mrt * uy * XKMPER, mrt * uz * XKMPER,
                           (mvt * ux + rvdot * vx) * VKMPERSEC,
                           (mvt * uy + rvdot * vy) * VKMPERSEC,
                           (mvt * uz + rvdot * vz) * VKMPERSEC);
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
