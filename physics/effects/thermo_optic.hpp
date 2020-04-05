/**
 * @file thermo_optic.hpp
 * @brief Thermo-optic effect in silicon photonic
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-03-11
 */

#ifndef OMNICLO_THERMO_OPTIC_HPP
#define OMNICLO_THERMO_OPTIC_HPP

#include "effects_base.hpp"

namespace effects {

    /**
     * Physical Realization for Thermo-Optic Effect
     */
    namespace thermo {

        /**
         * Constants for thermal-effect related materials
         */
        namespace constants {

            // Silicon
            const double THERMO_OPTIC_COEF_Si = 1.8e-4; //!< thermo-optic coefficient of Si [unit: \f$ K^{-1} \f$]
            const double THERMAL_CONDUCT_COEF_Si = 168.0;    //!< thermal conductivity coefficient of Si [unit: \f$ W \cdot K^{-1} m^{-1} \f$]
            const double THERMAL_EXPAND_COEF_Si = 2.5e-6;   //!< thermal expansion coefficient of Si [unit: \f$ K^{-1} \f$]

            // Silicon Dioxide
            const double THERMO_OPTIC_COEF_SiO2 = 0.1e-4; //!< thermo-optic coefficient of Si [unit: \f$ K^{-1} \f$]
            const double THERMAL_CONDUCT_COEF_SiO2 = 1.4;    //!< thermal conductivity coefficient of Si [unit: \f$ W \cdot K^{-1} m^{-1} \f$]
            const double THERMAL_EXPAND_COEF_SiO2 = 0.6e-6;   //!< thermal expansion coefficient of Si [unit: \f$ K^{-1} \f$]
        }

        /**
         * Calculate the change of refractive index of thermal-optic effect:
         * [ \f$ \Delta n = a\Delta T \f$ ]
         * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
         * @param TO_coef: thermal-optic coefficient [ \f$ a \f$ ]
         * @param delta_T: change of temperature in unit: *K* [ \f$ \Delta T \f$ ]
         * @return refractive index when thermo-optic effect occurs
         * @sa cal_thermo_optic_coef()
         */
        template <typename ValType>
        ValType cal_delta_n(ValType TO_coef, ValType delta_T) {
            static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                          std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

            return TO_coef * delta_T;
        }

        /**
         * Calculate the thermo-optic coefficient of thermo-optic effect:
         * [ \f$ a = \frac{\partial n}{\partial T} = \frac{(n^2 - 1)(n^2 + 2)}{6n} \cdot \left[ \frac{1}{\alpha} \frac{d\alpha}{dT} - 3\gamma \right] \f$ ]
         * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
         * @param n: refractive index at \f$ T_0 \f$ [ \f$ n_0 \f$ ]
         * @param avg_polar: average polarizability [ \f$ \alpha \f$ ]
         * @param partial_polar_T: differential average polarizability versus temperature [ \f$ \frac{\partial \alpha}{\partial T}  \f$ ]
         * @param linear_expand_coef: linear expansion coefficient [ \f$ \gamma \f$ ]
         * @return thermo-optic coefficient
         */
        template <typename ValType>
        ValType cal_thermo_optic_coef(ValType n, ValType avg_polar, ValType partial_polar_T, ValType linear_expand_coef) {
            static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                          std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

            ValType n2 = pow(n, 2);
            return (n2 - 1.0) * (n2 + 2.0) / (6.0 * n) * (1 / avg_polar * partial_polar_T - 3.0 * linear_expand_coef);
        }

        /**
         * Calculate the phase shift when light propagate in variant media:
         * [ \f$ \Delta \phi = \frac{2\pi}{\lambda_0} \cdot a L \Delta T \f$ ]
         * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
         * @param vacuum_wavelength: wavelength in unit: *micrometer* vacuum [ \f$ \lambda_0 \f$ ]
         * @param TO_coef: thermal-optic coefficient [ \f$ a \f$ ]
         * @param propagate_length: length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
         * @param delta_T: change of temperature in unit: *K* [ \f$ \Delta T \f$ ]
         * @return phase shift in radian
         * @sa cal_thermo_optic_coef()
         */
        template <typename ValType>
        ValType cal_phase_shift(ValType vacuum_wavelength, ValType TO_coef, ValType propagate_length, ValType delta_T) {
            static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                          std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

            return 2.0 * effects::constants::PI / vacuum_wavelength * TO_coef * propagate_length * delta_T;
        }

    }

}

#endif //OMNICLO_THERMO_OPTIC_HPP
