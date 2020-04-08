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
#include <utility>

#include "cells_base.hpp"

using Eigen::MatrixXd;

namespace cells {

    using effects::constants::i;

    /**
     * MicroRing Base Class Derived from CellBase
     *
     * Containing couple_matrix() and propagate_matrix() method
     */
    class MicroRingBase : public CellBase {
    public:
        MicroRingBase(size_t in_size, size_t out_size, std::string name): CellBase(in_size, out_size, std::move(name)) {}

    protected:
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
         * @param propagate_length12: propagation length of item \f$ H_{12} \f$
         * @param propagate_length21: propagation length of item \f$ H_{21} \f$
         * @return matrix \f$ H_{p} \f$
         */
        static Eigen::Matrix2cd propagate_matrix(const std::complex<double> &wave_vector,
                double propagate_length12, double propagate_length21) {
            Eigen::Matrix2cd m;
            std::complex<double> k = wave_vector;
            m(0, 0) = m(1, 1) = 0.0;
            m(0, 1) = exp(i * k * propagate_length12);
            m(1, 0) = exp(-i * k * propagate_length21);
            return m;
        }

    };

    /**
     * Serial-Coupled MicroRing derived from MicroRingBase.
     *
     * Available parameters: size_t(ring_number), vector<double>(radius, couple_length, self_couple_coef, cross_couple_coef),
     * vector<complex<double>>(n, delta_n).
     * @note couple_length, self_couple_coef, cross_couple_coef has (ring_number + 1) size because of topology!
     */
    class MicroRingSerial : public MicroRingBase {
    public:
        MicroRingSerial(const std::string &name, const std::map<std::string, std::any> &vars):
            MicroRingBase(2, 2, name) {
            device_vars = vars;
            auto ring_num = std::any_cast<size_t>(utils::map_setup_default(device_vars, "ring_number", size_t(1)));
            utils::map_setup_default(device_vars, "radius", std::vector<double>(ring_num, 6.506));
            utils::map_setup_default(device_vars, "couple_length", std::vector<double>(ring_num + 1, 3.5));
            utils::map_setup_default(device_vars, "self_couple_coef", std::vector<double>(ring_num + 1, 0.912));
            utils::map_setup_default(device_vars, "cross_couple_coef", std::vector<double>(ring_num + 1, 0.410));
            //utils::map_setup_default(device_vars, "n", std::vector<std::complex<double>>(ring_num, 2.5));
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
            //const auto &_n = std::any_cast<std::vector<std::complex<double>>>(device_vars["n"]);
            const auto &_delta_n = std::any_cast<std::vector<std::complex<double>>>(device_vars["delta_n"]);

            auto _k0 = 2.0 * PI * (effects::n_effective(wavelength) + _delta_n[0]) / wavelength + i * effects::propagating_loss();
            Matrix2cd m = couple_matrix(_self_couple[0], _cross_couple[0], _k0 ,_couple_length[0]);
            for(size_t cnt = 0; cnt < _ring_num; cnt++) {
                auto _n_eff = effects::n_effective(wavelength);
                auto _k = 2.0 * PI * (_n_eff + _delta_n[cnt]) / wavelength + i * effects::bending_loss(_radius[cnt]) + i * effects::propagating_loss();
                auto _prop_length = _radius[cnt] * PI;
                m = propagate_matrix(_k, _prop_length, _prop_length) * m;
                _k = 2.0 * PI * (_n_eff + _delta_n[cnt]) / wavelength + i * effects::propagating_loss();
                m = couple_matrix(_self_couple[cnt + 1], _cross_couple[cnt + 1], _k, _couple_length[cnt + 1]) * m;
            }
            Matrix2cd mt = transfer_matrix(m);
            E_out = mt * E_in;
        }

        /**
         * Coupler control strategy
         * @param ctrl_vars: control variables
         * @note Default ctrl_func provides thermo and electro control that pre-defined in effects
         * @note use delta_Nh, delta_Ne, delta_T in double to control the switches
         */
        void control(std::map<std::string, std::any> ctrl_vars) override {
            if(!ctrl_func) {
                if(utils::map_key_exist(ctrl_vars, "delta_Nh") && utils::map_key_exist(ctrl_vars, "delta_Ne")) {
                    const auto &_ring_num = std::any_cast<size_t>(device_vars["ring_number"]);
                    auto _delta_Nh = std::any_cast<double>(ctrl_vars["delta_Nh"]);
                    auto _delta_Ne = std::any_cast<double>(ctrl_vars["delta_Ne"]);
                    auto _delta_n = effects::electro::cal_delta_n<double>(_delta_Ne, _delta_Nh,
                            effects::electro::constants::COEF_P_1550, effects::electro::constants::COEF_Q_1550,
                            effects::electro::constants::COEF_R_1550, effects::electro::constants::COEF_S_1550);
                    device_vars["delta_n"] = std::vector<std::complex<double>>(_ring_num, _delta_n);
                } else if(utils::map_key_exist(ctrl_vars, "delta_T")) {
                    const auto &_ring_num = std::any_cast<size_t>(device_vars["ring_number"]);
                    auto _delta_T = std::any_cast<double>(ctrl_vars["delta_T"]);
                    auto _delta_n = effects::thermo::cal_delta_n<double>(effects::thermo::constants::THERMO_OPTIC_COEF_Si, _delta_T);
                    device_vars["delta_n"] = std::vector<std::complex<double>>(_ring_num, _delta_n);
                }
            } else {
                auto ctrl_targets = ctrl_func(ctrl_vars);
            }
        }

    protected:

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

        __deprecated double radius{};    //!< radius of MR-Mono
        __deprecated std::complex<double> kappa1, kappa2;    //!< coupling coefficients
        __deprecated std::complex<double> trans1, trans2;    //!< transmission coefficients
        __deprecated std::complex<double> refractive_index;  //!< refractive index of waveguide [ \f$ n \f$ ]
        __deprecated std::complex<double> delta_n;   //!< change of refractive index of waveguide [ \f$ \Delta n \f$ ]

    };



    /**
     * MicroRing Resonator with Cross Waveguide
     */
    class MicroRingCross : public MicroRingBase {
    public:
        MicroRingCross(const std::string &name, const std::map<std::string, std::any> &vars):
                MicroRingBase(2, 2, name) {
            device_vars = vars;
            auto ring_num = std::any_cast<size_t>(utils::map_setup_default(device_vars, "ring_number", size_t(1)));
            utils::map_setup_default(device_vars, "radius", std::vector<double>(ring_num, 6.502));
            utils::map_setup_default(device_vars, "couple_length", std::vector<double>(ring_num + 1, 3.5));
            utils::map_setup_default(device_vars, "self_couple_coef", std::vector<double>(ring_num + 1, 0.912));
            utils::map_setup_default(device_vars, "cross_couple_coef", std::vector<double>(ring_num + 1, 0.410));
            //utils::map_setup_default(device_vars, "n", std::vector<std::complex<double>>(ring_num, 1.0));
            utils::map_setup_default(device_vars, "delta_n", std::vector<std::complex<double>>(ring_num, 0.0));
            utils::map_setup_default(device_vars, "insertion_loss_dB", double(-0.11));
            utils::map_setup_default(device_vars, "crosstalk_dB", double(-45.2));
        }
        ~MicroRingCross() override = default;

        void update(double wavelength) override {
            using Eigen::Matrix2cd;
            using effects::constants::PI;
            const auto &_ring_num = std::any_cast<size_t>(device_vars["ring_number"]);
            const auto &_radius = std::any_cast<std::vector<double>>(device_vars["radius"]);
            const auto &_couple_length = std::any_cast<std::vector<double>>(device_vars["couple_length"]);
            const auto &_self_couple = std::any_cast<std::vector<double>>(device_vars["self_couple_coef"]);
            const auto &_cross_couple = std::any_cast<std::vector<double>>(device_vars["cross_couple_coef"]);
            //const auto &_n = std::any_cast<std::vector<std::complex<double>>>(device_vars["n"]);
            const auto &_delta_n = std::any_cast<std::vector<std::complex<double>>>(device_vars["delta_n"]);
            auto insertion_loss_dB = std::any_cast<double>(device_vars["insertion_loss_dB"]);
            auto crosstalk_dB = std::any_cast<double>(device_vars["crosstalk_dB"]);

            auto _k0 = 2.0 * PI * (effects::n_effective(wavelength) + _delta_n[0]) / wavelength + i * effects::propagating_loss();
            Matrix2cd m_ring = couple_matrix(_self_couple[0], _cross_couple[0], _k0 ,_couple_length[0]);
            for(size_t cnt = 0; cnt < _ring_num; cnt++) {
                auto _n_eff = effects::n_effective(wavelength);
                auto _k = 2.0 * PI * (_n_eff + _delta_n[cnt]) / wavelength + i * effects::bending_loss(_radius[cnt]) + i * effects::propagating_loss();
                auto _prop_length = _radius[cnt] * PI;
                m_ring = propagate_matrix(_k, _prop_length, _prop_length) * m_ring;
                _k = 2.0 * PI * (_n_eff + _delta_n[cnt]) / wavelength + i * effects::propagating_loss();
                m_ring = couple_matrix(_self_couple[cnt + 1], _cross_couple[cnt + 1], _k, _couple_length[cnt + 1]) * m_ring;
            }

            Matrix2cd m_cross;
            m_cross(0, 0) = m_cross(1, 1) = exp(insertion_loss_dB / 10.0);
            m_cross(1, 0) = m_cross(0, 1) = exp(crosstalk_dB / 10.0);
            Matrix2cd mt = transfer_matrix(m_ring, m_cross);
            E_out = mt * E_in;
        }

        void control(std::map<std::string, std::any> ctrl_vars) override {

        }

    protected:
        static Eigen::Matrix2cd transfer_matrix(const Eigen::Matrix2cd &m_ring, const Eigen::Matrix2cd &m_cross) {
            Eigen::Matrix2cd m;
            auto denom = m_ring(0, 1) - m_cross(1, 0);
            m(0, 0) = -m_cross(1, 0) * m_ring(0, 0) / denom;
            m(0, 1) = m_cross(1, 1) * m_ring(0, 1) / denom;
            m(1, 0) =
                    (m_ring(0, 1) * m_ring(1, 0) - m_ring(1, 0) * m_cross(1, 0) - m_ring(0, 0) * m_ring(1, 1)) / denom;
            m(1, 1) = m_cross(1, 1) * m_ring(1, 1) / denom;
            return m;
        }
    };

}

#endif //OMNICLO_MICRORING_HPP
