/**
 * @file electro_optic.hpp
 * @brief Electro-optic effect in silicon photonic
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-03-11
 */

#ifndef OMNICLO_ELECTRO_OPTIC_HPP
#define OMNICLO_ELECTRO_OPTIC_HPP

#include "effects_base.hpp"

namespace effects {

    /**
     * Physical Realization for Electro-Optic Effect
     */
    namespace electro {

         /**
          * Constants for electro-effect related materials
          */
         namespace constants {

             // absorption and refractive index coefficients in wavelength of 1310 nm
             const double COEF_A_1310 = 3.48e-21;
             const double COEF_B_1310 = 1.229;
             const double COEF_C_1310 = 1.02e-19;
             const double COEF_D_1310 = 1.089;
             const double COEF_P_1310 = 2.98e-22;
             const double COEF_Q_1310 = 1.016;
             const double COEF_R_1310 = 1.25e-18;
             const double COEF_S_1310 = 0.835;

             // absorption and refractive index coefficients in wavelength of 1550 nm
             const double COEF_A_1550 = 8.88e-21;
             const double COEF_B_1550 = 1.167;
             const double COEF_C_1550 = 5.84e-20;
             const double COEF_D_1550 = 1.109;
             const double COEF_P_1550 = 5.40e-22;
             const double COEF_Q_1550 = 1.011;
             const double COEF_R_1550 = 1.53e-18;
             const double COEF_S_1550 = 0.838;
         }

         /**
          * Calculate the change of refractive index of electro-optic effect:
          * [ \f$ \Delta n = -\frac{n^3 \gamma }{2} E_0 \f$ ]
          * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
          * @param EO_coef: electro-optic coefficient [ \f$ \gamma \f$ ]
          * @param n: permittivity of the materials [ \f$ n \f$ ]
          * @param E_field: electric field takes effects on materials [ \f$ E_0 \f$ ]
          * @return refractive index when electro-optic effect occurs
          */
         template <typename ValType>
         ValType cal_delta_n(ValType EO_coef, ValType n, ValType E_field) {
             static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                          std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

             return -(pow(n, 3) * EO_coef * E_field) / 2;
         }

        /**
         * Calculate the change of refractive index of electro-optic effect:
         * [ \f$ \Delta n = -[p * \Delta N_e^{q} + r * \Delta N_h^s] \f$ ]
         * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
         * @param delta_Ne: change of the concentration of electrons in unit: \f$ cm^{-3} \f$ [ \f$ \Delta N_e \f$ ]
         * @param delta_Nh: change of the concentration of holes in unit: \f$ cm^{-3} \f$ [ \f$ \Delta N_h \f$ ]
         * @param coef_p: coefficient p [ \f$ p(\lambda) \f$ ]
         * @param coef_q: coefficient q [ \f$ q(\lambda) \f$ ]
         * @param coef_r: coefficient r [ \f$ r(\lambda) \f$ ]
         * @param coef_s: coefficient s [ \f$ s(\lambda) \f$ ]
         * @return refractive index when electro-optic effect occurs
         */
        template <typename ValType>
        ValType cal_delta_n(ValType delta_Ne, ValType delta_Nh,
                ValType coef_p, ValType coef_q, ValType coef_r, ValType coef_s) {
            static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                          std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

            return -(coef_p * pow(delta_Ne, coef_q) + coef_r * pow(delta_Nh, coef_s));
        }

        /**
         * Calculate the change of absorption loss coefficient of electro-optic effect:
         * [ \f$ \Delta \alpha = [a * \Delta N_e^{b} + c * \Delta N_h^d] \f$ ]
         * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
         * @param delta_Ne: change of the concentration of electrons in unit: \f$ cm^{-3} \f$ [ \f$ \Delta N_e \f$ ]
         * @param delta_Nh: change of the concentration of holes in unit: \f$ cm^{-3} \f$ [ \f$ \Delta N_h \f$ ]
         * @param coef_a: coefficient a [ \f$ a(\lambda) \f$ ]
         * @param coef_b: coefficient b [ \f$ b(\lambda) \f$ ]
         * @param coef_c: coefficient c [ \f$ c(\lambda) \f$ ]
         * @param coef_d: coefficient d [ \f$ d(\lambda) \f$ ]
         * @return refractive index when electro-optic effect occurs
         */
        template <typename ValType>
        ValType cal_delta_alpha(ValType delta_Ne, ValType delta_Nh,
                ValType coef_a, ValType coef_b, ValType coef_c, ValType coef_d) {
            static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                          std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

            return (coef_a * pow(delta_Ne, coef_b) + coef_c * pow(delta_Nh, coef_d));
        }

         /**
          * Calculate the phase shift when light propagate in variant media:
          * [ \f$ \Delta \phi = -\frac{\pi n^3 \gamma L}{\lambda_0} E_0 \f$ ]
          * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
          * @param vacuum_wavelength: wavelength in unit: *micrometer* vacuum [ \f$ \lambda_0 \f$ ]
          * @param EO_coef: electro-optic coefficient [ \f$ \gamma \f$ ]
          * @param propagate_length: length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
          * @param n: permittivity of the materials [ \f$ n \f$ ]
          * @param E_field: electric field takes effects on materials [ \f$ E_0 \f$ ]
          * @return phase shift in radian
          */
         template <typename ValType>
         ValType cal_phase_shift(ValType vacuum_wavelength, ValType EO_coef, ValType propagate_length, ValType n, ValType E_field) {
             static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                          std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

             return -(effects::constants::PI * pow(n, 3) * EO_coef * propagate_length * E_field) / vacuum_wavelength;
         }

         /**
          * Calculate the nonlinear coefficients:
          * [ \f$ \gamma = \frac{k_0 n_2}{S} = \frac{2\pi n_2}{\lambda_0 S} \f$ ]
          * @tparam ValType: Value Type (float, double, std::valarray<float>, std::valarray<double>)
          * @param vacuum_wavelength: wavelength in unit: *micrometer* vacuum [ \f$ \lambda_0 \f$ ]
          * @param n2: second order of nonlinear index of refractive [ \f$ n_2 \f$ ]
          * @param waveguide_sectional_area: sectional area of waveguide in unit: \f$ \mu m^2 \f$ [ \f$ S \f$ ]
          * @return nonlinear coefficients
          */
         template <typename ValType>
         ValType cal_nonlinear_coef(ValType vacuum_wavelength, ValType n2, ValType waveguide_sectional_area) {
             static_assert(std::is_same<ValType, double>::value || std::is_same<ValType, float>::value ||
                           std::is_same<ValType, std::valarray<double>>::value || std::is_same<ValType, std::valarray<float>>::value);

             return (2 * effects::constants::PI * n2) / (vacuum_wavelength * waveguide_sectional_area);
         }
    }
}

#endif //OMNICLO_ELECTRO_OPTIC_HPP
