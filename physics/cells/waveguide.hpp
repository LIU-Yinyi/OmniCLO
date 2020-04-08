/**
 * @file waveguide.hpp
 * @brief Relization of Waveguide Generalized Cell
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_WAVEGUIDE_HPP
#define OMNICLO_WAVEGUIDE_HPP

// Custom Library Header
#include "cells_base.hpp"

namespace cells {

    using effects::constants::i;

    /**
     * Waveguide derived from CellBase.
     */
    class Waveguide : public CellBase {
    public:
        /**
         * Constructor with Waveguide's properties
         * @param name: name of device cell
         * @param vars: parameters of device cell
         */
        Waveguide(const std::string &name, const std::map<std::string, std::any> &vars):
            CellBase(1, 1, name) {
            device_vars = vars;
            //utils::map_setup_default(device_vars, "n", std::complex<double>(1, 0));
            utils::map_setup_default(device_vars, "delta_n", std::complex<double>(0, 0));
            utils::map_setup_default(device_vars, "propagate_length", double(10.0));
        }
        ~Waveguide() override = default;

        /**
         * Waveguide update strategy
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         */
        void update(double wavelength) override {
            using effects::constants::PI;
            //const auto &_refractive_index = std::any_cast<std::complex<double>>(device_vars["n"]);
            auto _refractive_index = effects::n_effective(wavelength);
            const auto &_delta_n = std::any_cast<std::complex<double>>(device_vars["delta_n"]);
            const auto &_propagate_length = std::any_cast<double>(device_vars["propagate_length"]);
            std::complex<double> k = 2.0 * PI * (_refractive_index + _delta_n) / wavelength + i * effects::propagating_loss();
            E_out = exp(i * k * _propagate_length) * E_in;
        }

        /**
         * Waveguide control strategy
         * @param ctrl_vars: control variables
         * @note CtrlFunc return value needs 1 element for delta_n [ \f$ \Delta n \f$ ]
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
        __deprecated std::complex<double> refractive_index;  //!< refractive index of waveguide [ \f$ n \f$ ]
        __deprecated std::complex<double> delta_n;   //!< change of refractive index of waveguide [ \f$ \Delta n \f$ ]
        __deprecated double propagate_length{};  //!< length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
    };



    /**
     * Waveguide derived from CellBase.
     */
    class WaveguideCross : public CellBase {
    public:
        /**
         * Constructor with Waveguide Cross's properties
         * @param name: name of device cell
         * @param vars: parameters of device cell
         * @note cite: Celo, D., et al. "Low-loss waveguide crossings for photonic integrated circuits on SOI technology."
         * 11th International Conference on Group IV Photonics (GFP). IEEE, 2014.
         */
        WaveguideCross(const std::string &name, const std::map<std::string, std::any> &vars) :
                CellBase(2, 2, name) {
            device_vars = vars;
            utils::map_setup_default(device_vars, "insertion_loss_dB", double(-0.11));
            utils::map_setup_default(device_vars, "crosstalk_dB", double(-45.2));
        }

        ~WaveguideCross() override = default;

        /**
         * Waveguide update strategy
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         */
        void update(double wavelength) override {
            auto insertion_loss_dB = std::any_cast<double>(device_vars["insertion_loss_dB"]);
            auto crosstalk_dB = std::any_cast<double>(device_vars["crosstalk_dB"]);
            Eigen::Matrix2cd m;
            double Li = exp(insertion_loss_dB / 10.0);
            double Xt = exp(crosstalk_dB / 10.0);
            m(0, 0) = m(1, 1) = Li;
            m(1, 0) = m(0, 1) = Xt;
            E_out = m * E_in;
        }

        /**
         * Waveguide control strategy
         * @param ctrl_vars: control variables
         * @note CtrlFunc return value needs 1 element for delta_n [ \f$ \Delta n \f$ ]
         */
        void control(std::map<std::string, std::any> ctrl_vars) override {
            
        }
    };
}

#endif //OMNICLO_WAVEGUIDE_HPP
