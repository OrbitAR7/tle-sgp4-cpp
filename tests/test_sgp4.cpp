// Checks the C++ SGP4 against reference states from the `sgp4` Python package
// (Vallado's implementation). Regenerate with scripts/make_reference.py.
#include "sgp4_propagator.hpp"
#include <cmath>
#include <cstdio>

struct Ref {
    const char* name; const char* l1; const char* l2;
    double t, x, y, z, vx, vy, vz;
};

static const Ref REFS[] = {
#include "reference_states.inc"
};

static int failures = 0;
static void check(bool ok, const char* what) {
    std::printf("%s  %s\n", ok ? "pass" : "FAIL", what);
    if (!ok) ++failures;
}

int main() {
    double worst_r = 0.0, worst_v = 0.0;
    for (const Ref& r : REFS) {
        sgp4::SGP4 sat(r.l1, r.l2, r.name);
        sgp4::StateVector s = sat.propagate(r.t);
        double dr = std::sqrt(std::pow(s.x - r.x, 2) + std::pow(s.y - r.y, 2) + std::pow(s.z - r.z, 2));
        double dv = std::sqrt(std::pow(s.vx - r.vx, 2) + std::pow(s.vy - r.vy, 2) + std::pow(s.vz - r.vz, 2));
        worst_r = std::fmax(worst_r, dr);
        worst_v = std::fmax(worst_v, dv);
    }
    std::printf("max position difference vs reference: %.3e m\n", worst_r * 1e3);
    std::printf("max velocity difference vs reference: %.3e mm/s\n", worst_v * 1e6);
    check(worst_r < 1e-6, "position within 1 mm of reference SGP4");
    check(worst_v < 1e-9, "velocity within 1 um/s of reference SGP4");

    // Epoch: 2024 day 1.5 is 2024-01-01 12:00 UTC = JD 2460311.0
    sgp4::TLE tle = sgp4::parseTLE(REFS[0].l1, REFS[0].l2);
    check(std::fabs(tle.epoch_jd - 2460311.0) < 1e-9, "TLE epoch Julian date");

    // Deep-space objects are rejected rather than propagated with the wrong model
    sgp4::SGP4 gps("1 24876U 97035A   24001.50000000 -.00000020  00000-0  00000+0 0  9999",
                   "2 24876  55.7000  50.0000 0100000 270.0000  90.0000  2.00560000100000");
    bool threw = false;
    try { gps.propagate(0.0); } catch (const std::domain_error&) { threw = true; }
    check(gps.isDeepSpace() && threw, "deep-space TLE detected and rejected");

    std::printf("%s\n", failures ? "SOME TESTS FAILED" : "all tests passed");
    return failures ? 1 : 0;
}
