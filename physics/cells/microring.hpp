/**
 * @file microring.hpp
 * @brief Relization of Micro-Ring Resonator Switching Cell
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_MICRORING_HPP
#define OMNICLO_MICRORING_HPP

// Custom Library Header
#include "cells_base.hpp"

using Eigen::MatrixXd;

namespace cells {

    using effects::constants::i;

    /**
     * Serial-Coupled MicroRing derived from CellBase.
     *
     * Available parameters: size_t(ring_number), vector<double>(radius, couple_length, self_couple_coef, cross_couple_coef),
     * vector<complex<double>>(n, delta_n).
     * @note couple_length, self_couple_coef, cross_couple_coef has (ring_number + 1) size because of topology!
     */
    class MicroRingSerial : public CellBase {
    public:
        MicroRingSerial(const std::string &name, const std::map<std::string, std::any> &vars):
            CellBase(2, 2, name) {
            device_vars = vars;
            auto ring_num = std::any_cast<size_t>(utils::map_setup_default(device_vars, "ring_number", size_t(1)));
            utils::map_setup_default(device_vars, "radius", std::vector<double>(ring_num, 6.5));
            utils::map_setup_default(device_vars, "couple_length", std::vector<double>(ring_num + 1, 3.5));
            utils::map_setup_default(device_vars, "self_couple_coef", std::vector<double>(ring_num + 1, 0.912));
            utils::map_setup_default(device_vars, "cross_couple_coef", std::vector<double>(ring_num + 1, 0.410));
            utils::map_setup_default(device_vars, "n", std::vector<std::complex<double>>(ring_num, 1.0));
            utils::map_setup_default(device_vars, "delta_n", std::vector<std::complex<double>>(ring_num, 0.0));
        }
        ~MicroRingSerial() override = default;

        /**
         * Coupler update strategy
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         */
        void update(double wavelength) override {
            using Eigen::Matrix2cd;
            using effects::constants::PI;
            const auto &_ring_num = std::any_cast<size_t>(device_vars["ring_number"]);
            const auto &_radius = std::any_cast<std::vector<double>>(device_vars["radius"]);
            const auto &_couple_length = std::any_cast<std::vector<double>>(device_vars["couple_length"]);
            const auto &_self_couple = std::any_cast<std::vector<double>>(device_vars["self_couple_coef"]);
            const auto &_cross_couple = std::any_cast<std::vector<double>>(device_vars["cross_couple_coef"]);
            const auto &_n = std::any_cast<std::vector<std::complex<double>>>(device_vars["n"]);
            const auto &_delta_n = std::any_cast<std::vector<std::complex<double>>>(device_vars["delta_n"]);
            Matrix2cd m = couple_matrix(_self_couple[0], _cross_couple[0],
                    2.0 * PI * (_n[0] + _delta_n[0]) / wavelength, _couple_length[0]);
            for(size_t cnt = 0; cnt < _ring_num; cnt++) {
                auto _k = 2.0 * PI * (_n[cnt] + _delta_n[cnt]) / wavelength;
                m = propagate_matrix(_k, _radius[cnt] * PI) * m;
                m = couple_matrix(_self_couple[cnt + 1], _cross_couple[cnt + 1], _k, _couple_length[cnt + 1]) * m;
            }
            Matrix2cd mt = transfer_matrix(m);
            E_out = mt * E_in;
        }

        /**
         * Coupler control strategy
         * @param ctrl_vars: control variables
         * @note CtrlFunc return value needs 1 element for delta_beta [ \f$ \Delta \beta \f$ ]
         */
        void control(std::map<std::string, std::any> ctrl_vars) override {
            if(!ctrl_func) {

            } else {
                auto ctrl_targets = ctrl_func(ctrl_vars);
            }
        }

    protected:
        static double bending_loss(double ring_radius) {
            double _alpha = 0.0;
            if(ring_radius <= 20.0) {
                _alpha = 0.03124 * pow(ring_radius, -3.02056);
            }
            return _alpha;
        }

        /**
         * Acquire the Couple Matrix: \n
         * \f$ \bmatrix{ E_a \\ E_d} = H_{c} \bmatrix{ E_i \\ E_t } \f$ \n
         * @param self_couple_coef: self-coupled coefficient
         * @param cross_couple_coef: cross-coupled coefficient
         * @param wave_vector: wave vector \f$ k = \beta + i\alpha = \frac{2\pi}{\lambda}(n + \Delta n) \f$
         * @param couple_length: coupled length of nearly straight section
         * @return matrix \f$ H_{c} \f$
         */
        static Eigen::Matrix2cd couple_matrix(double self_couple_coef, double cross_couple_coef,
                const std::complex<double> &wave_vector, double couple_length) {
            Eigen::Matrix2cd m;
            double r = self_couple_coef;
            double t = cross_couple_coef;
            std::complex<double> k = wave_vector;
            m(0, 0) = i * r / t;
            m(0, 1) = - i * exp(-i * k * couple_length) / t;
            m(1, 0) = i * (t * t + r * r) / t * exp(i * k * couple_length);
            m(1, 1) = -i * r / t;
            return m;
        }

        /**
         * Acquire the Propagate Matrix: \n
         * \f$ H_{p} = \bmatrix{ 0 & e^{ikL} \\ e^{-ikL} & 0 } \f$
         * @param wave_vector
         * @param propagate_length
         * @return matrix \f$ H_{p} \f$
         */
        static Eigen::Matrix2cd propagate_matrix(const std::complex<double> &wave_vector, double propagate_length) {
            Eigen::Matrix2cd m;
            std::complex<double> k = wave_vector;
            m(0, 0) = m(1, 1) = 0.0;
            m(0, 1) = exp(i * k * propagate_length);
            m(1, 0) = exp(-i * k * propagate_length);
            return m;
        }

        /**
         * Acquire transfer matrix from TYPE-1 to TYPE-2 \n
         * \f$ \bmatrix{ E_a \\ E_d} = H \bmatrix{ E_i \\ E_t } \f$ \n
         * \f$ \bmatrix{ E_t \\ E_d} = H' \bmatrix{ E_i \\ E_a } \f$ \n
         * @param ori_m: origin transfer matrix \f$ H \f$
         * @return transformed transfer matrix \f$ H' \f$
         */
        static Eigen::Matrix2cd transfer_matrix(const Eigen::Matrix2cd &ori_m) {
            Eigen::Matrix2cd m;
            m(0, 0) = -ori_m(0, 0) / ori_m(0, 1);
            m(0, 1) = 1.0 / ori_m(0, 1);
            m(1, 0) = (ori_m(0, 1) * ori_m(1, 0) - ori_m(1, 1) * ori_m(0, 0)) / ori_m(0, 1);
            m(1, 1) = ori_m(1, 1) / ori_m(0, 1);
            return m;
        }

    protected:
        __deprecated double radius{};    //!< radius of MR-Mono
        __deprecated std::complex<double> kappa1, kappa2;    //!< coupling coefficients
        __deprecated std::complex<double> trans1, trans2;    //!< transmission coefficients
        __deprecated std::complex<double> refractive_index;  //!< refractive index of waveguide [ \f$ n \f$ ]
        __deprecated std::complex<double> delta_n;   //!< change of refractive index of waveguide [ \f$ \Delta n \f$ ]

    };



    class MicroRingCross : public CellBase {
    public:
        MicroRingCross(const std::string &name, const std::map<std::string, std::any> &vars):
                CellBase(2, 2, name) {
            device_vars = vars;
            auto ring_num = std::any_cast<size_t>(utils::map_setup_default(device_vars, "ring_number", size_t(1)));
            utils::map_setup_default(device_vars, "radius", std::vector<double>(ring_num, 6.5));
            utils::map_setup_default(device_vars, "couple_length", std::vector<double>(ring_num + 1, 3.5));
            utils::map_setup_default(device_vars, "self_couple_coef", std::vector<double>(ring_num + 1, 0.912));
            utils::map_setup_default(device_vars, "cross_couple_coef", std::vector<double>(ring_num + 1, 0.410));
            utils::map_setup_default(device_vars, "n", std::vector<std::complex<double>>(ring_num, 1.0));
            utils::map_setup_default(device_vars, "delta_n", std::vector<std::complex<double>>(ring_num, 0.0));
        }
        ~MicroRingCross() override = default;

        void update(double wavelength) override {

        }

        void control(std::map<std::string, std::any> ctrl_vars) override {

        }
    };

}
#endif //OMNICLO_MICRORING_HPP
