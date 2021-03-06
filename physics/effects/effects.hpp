/**
 * @file effects.hpp
 * @brief Integral effects in silicon photonic
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-03-11
 */

#ifndef OMNICLO_EFFECTS_HPP
#define OMNICLO_EFFECTS_HPP

#include "thermo_optic.hpp"
#include "electro_optic.hpp"

namespace effects {

    /**
     * Get the propagation loss of waveguide
     * @param func_ptr: customize function pointer used for experiment data manager
     * @return absorption coefficients \f$ alpha \f$
     */
    inline double propagating_loss(double (*func_ptr)() = nullptr) { return (func_ptr != nullptr) ? func_ptr() : 3.0e-5; }

    /**
     * Get the bending loss of ring structure
     * @param ring_radius: ring radius in unit: *micrometer*
     * @param func_ptr: customize function pointer used for experiment data manager
     * @return absorption coefficients \f$ alpha \f$
     */
    double bending_loss(double ring_radius, double(*func_ptr)(double) = nullptr) {
        double _alpha = 0.0;
        if(!func_ptr) {
            if (ring_radius <= 20.0 && ring_radius >= 2.0) {
                _alpha = 0.03124 * pow(ring_radius, -3.02056);
            }
        } else {
            _alpha = func_ptr(ring_radius);
        }
        return _alpha;
    }

    /**
     * Get the effective refractive index
     * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
     * @param func_ptr: customize function pointer used for experiment data manager
     * @return effective refractive index
     */
    std::complex<double> n_effective(double wavelength, std::complex<double>(*func_ptr)(double) = nullptr) {
        using effects::constants::PI;
        std::complex<double> _retval;
        if(!func_ptr) {
            _retval = -0.9885 * wavelength + 4.0347;
        } else {
            _retval = func_ptr(wavelength);
        }
        return _retval;
    }

}

#endif //OMNICLO_EFFECTS_HPP
