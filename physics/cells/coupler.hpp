/**
 * @file coupler.hpp
 * @brief Relization of Coupler Cell
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_COUPLER_HPP
#define OMNICLO_COUPLER_HPP

// Custom Library Header
#include "cells_base.hpp"

namespace cells {

    using effects::constants::i;

    /**
     * Coupler derived from CellBase.
     *
     * When \f$ \kappa_12 = \kappa_21 \f$ means Symmetrical Coupler, \n
     * When \f$ \kappa_12 \neq \kappa_21 \f$ means Asymmetrical Coupler.
     */
    class Coupler : public CellBase {
    public:
        /**
         * Constructor with Coupler's properties
         * @param name: name of device cell
         * @param vars: parameters of device cell
         */
        Coupler(const std::string &name, const std::map<std::string, std::any> &vars):
            CellBase(2, 2, name) {
            utils::map_setup_default(device_vars, "n", std::complex<double>(1, 0));
            utils::map_setup_default(device_vars, "delta_n", std::complex<double>(0, 0));
            utils::map_setup_default(device_vars, "propagate_length", double(10.0));
            utils::map_setup_default(device_vars, "kappa12", double(1.0));
            utils::map_setup_default(device_vars, "kappa21", double(1.0));
        }
        ~Coupler() override = default;

        /**
         * Coupler update strategy
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         */
        void update(double wavelength) override {
            using Eigen::Matrix2cd;
            using effects::constants::PI;
            const auto &_kappa12 = std::any_cast<double>(device_vars["kappa12"]);
            const auto &_kappa21 = std::any_cast<double>(device_vars["kappa21"]);
            //const auto &_refractive_index = std::any_cast<std::complex<double>>(device_vars["n"]);
            const auto &_delta_n = std::any_cast<std::complex<double>>(device_vars["delta_n"]);
            const auto &_propagate_length = std::any_cast<double>(device_vars["propagate_length"]);
            Matrix2cd m;
            double delta_beta = 2.0 * PI * _delta_n.real() / wavelength;
            double kappa = sqrt(_kappa12 * _kappa21);
            double K = sqrt(pow(delta_beta/2.0, 2.0) + pow(kappa, 2.0));
            double Kz = K * _propagate_length;
            double cosKz = cos(Kz);
            double sinKz = sin(Kz);
            std::complex<double> exp_ipz = exp(i * delta_beta * _propagate_length / 2.0);
            std::complex<double> exp_ipzn = exp(-i * delta_beta * _propagate_length / 2.0);
            std::complex<double> d1 = exp_ipz * (cosKz - i * delta_beta * sinKz / (2 * K));
            std::complex<double> d2 = i * exp_ipzn * (kappa / K) * sinKz;
            m(0, 0) = d1;
            m(1, 1) = std::conj(d1);
            m(0, 1) = m(1, 0) = d2;
            E_out = m * E_in;
        }

        /**
         * Coupler control strategy
         * @param ctrl_vars: control variables
         * @note CtrlFunc return value needs 1 element for delta_beta [ \f$ \Delta \beta \f$ ]
         */
        void control(std::map<std::string, std::any> ctrl_vars) override {
            auto &_delta_n = std::any_cast<std::complex<double>&>(device_vars["delta_n"]);
            if(!ctrl_func) {
                if(utils::map_key_exist(ctrl_vars, "delta_n")) { _delta_n = std::any_cast<std::complex<double>>(ctrl_vars["delta_n"]); }
                else { _delta_n = 0.0; }
            } else {
                //TODO: update control function and variables
                auto ctrl_targets = ctrl_func(ctrl_vars);
                if(ctrl_targets.size() == 1) { _delta_n = std::any_cast<std::complex<double>>(ctrl_targets[0]); }
                else { _delta_n = 0.0; }
            }
        }

    protected:
        __deprecated double kappa1{}, kappa2{};  //!< coupling coefficients of two line waveguides
        __deprecated std::complex<double> refractive_index;  //!< refractive index of waveguide [ \f$ n \f$ ]
        __deprecated std::complex<double> delta_n;   //!< change of refractive index of waveguide [ \f$ \Delta n \f$ ]
        __deprecated double propagate_length{};    //!< length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
    };
}

#endif //OMNICLO_COUPLER_HPP
