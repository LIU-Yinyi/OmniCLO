/**
 * @file effects_base.hpp
 * @brief Base of physical effects in silicon photonic
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-03-11
 */

#ifndef OMNICLO_EFFECTS_BASE_HPP
#define OMNICLO_EFFECTS_BASE_HPP

// C++ Built-in Header
#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <random>
#include <utility>
#include <initializer_list>

// Scientific Library Header
#include <eigen3/Eigen/Dense>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

namespace effects {

    namespace constants {
        constexpr double PI = 3.1415926535897932384626433827950288; //!< ratio of a circle's circumference to its diameter
        constexpr std::complex<double> i = std::complex<double>(0.0, 1.0);  //!< imaginary part of complex number
    }

    /**
     * Calculate the refractive index from beta and alpha
     * @param beta: propagation constants [ \f$ \beta \f$ ]
     * @param alpha: absorption constants [ \f$ \alpha \f$ ]
     * @return refractive index in complex [ \f$ n \f$ ]
     */
    std::complex<double> cal_n(double beta, double alpha) {
        return std::complex<double>(beta, alpha);
    }

    /**
      * Calculate the change of propagation constant:
      * [ \f$ \Delta \beta = \frac{2\pi \Delta n}{\lambda_0} \f$ ]
      * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
      * @param delta_n: change of refractive index [ \f$ \Delta n \f$ ]
      * @param vacuum_wavelength: wavelength in unit: *micrometer* vacuum [ \f$ \lambda_0 \f$ ]
      * @return change of propagation constant
      * @sa cal_delta_n()
      */
    template <typename ValType>
    ValType cal_delta_beta(ValType delta_n, ValType vacuum_wavelength) {
        static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                      std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

        return (2.0 * effects::constants::PI * delta_n) / vacuum_wavelength;
    }

    /**
      * Calculate the phase shift when light propagate in variant media:
      * [ \f$ \Delta \phi = L \cdot \Delta\beta \f$ ]
      * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
      * @param delta_beta: the change of propagation constant [ \f$ \Delta \beta \f$ ]
      * @param propagate_length: length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
      * @return phase shift in radian
      */
    template <typename ValType>
    ValType cal_phase_shift(ValType delta_beta, ValType propagate_length) {
        static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                      std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

        return delta_beta * propagate_length;
    }

}

#endif //OMNICLO_EFFECTS_BASE_HPP
