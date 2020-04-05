/**
 * @file microring.hpp
 * @brief Relization of Micro-Ring Resonator Switching Cell
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-04-03
 */

#ifndef OMNICLO_MICRORING_HPP
#define OMNICLO_MICRORING_HPP

// Custom Library Header
#include "cells_base.hpp"

using Eigen::MatrixXd;

namespace cells {

    using effects::constants::i;

    class MicroRingArray : public CellBase {
    public:
        MicroRingArray(std::complex<double> n, double radius,
                double kappa1, double kappa2, double trans1, double trans2): CellBase(2, 2),
            refractive_index(n), delta_n(0.0), radius(radius),
            kappa1(kappa1), kappa2(kappa2), trans1(trans1), trans2(trans2) {}
        ~MicroRingArray() override = default;

        /**
         * Print the model properties and the protected or private parameters of the cell.
         */
        void print() override {

        }

        void update(double wavelength) override {

        }

        void control(std::vector<double> ctrl_vars) override {

        }

        /**
         * Update the parameters of cell properties that availbale for optimization.
         * @param optim_map: former key for availableparameter names, latter key for value
         * @note a complex number or higher dimension should be divided into several double
         * @note available parameters:
         */
        void optimize(const std::map<std::string, double> &optim_map) override {
            for(const auto &item : optim_map) {
                if(item.first == "n_real") { refractive_index.real(item.second); }
                else if(item.first == "n_imag") { refractive_index.imag(item.second); }
                else if(item.first == "delta_n_real") { delta_n.real(item.second); }
                else if(item.first == "delta_n_imag") { delta_n.imag(item.second); }
                else if(item.first == "radius") { radius = item.second; }
                else if(item.first == "kappa1") { kappa1 = item.second; }
                else if(item.first == "kappa2") { kappa2 = item.second; }
                else if(item.first == "trans1") { trans1 = item.second; }
                else if(item.first == "trans2") { trans2 = item.second; }
                else { std::cout << YELLOW << "[WARN] {cells->MicroRing-Mono/optimize} undefined optimization parameter: <"
                    << CYAN << item.first << YELLOW << ">." << RESET << std::endl; }
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

    protected:
        double radius;    //!< radius of MR-Mono
        std::complex<double> kappa1, kappa2;    //!< coupling coefficients
        std::complex<double> trans1, trans2;    //!< transmission coefficients
        std::complex<double> refractive_index;  //!< refractive index of waveguide [ \f$ n \f$ ]
        std::complex<double> delta_n;   //!< change of refractive index of waveguide [ \f$ \Delta n \f$ ]

    };

}
#endif //OMNICLO_MICRORING_HPP
