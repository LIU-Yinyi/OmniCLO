/**
 * @file terminator.hpp
 * @brief Relization of Terminator Generalized Cell
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_TERMINATOR_HPP
#define OMNICLO_TERMINATOR_HPP

// Custom Library Header
#include "cells_base.hpp"

namespace cells {

    /**
     * Input Port Derived from CellBase
     */
    class InputPort: public CellBase {
    public:
        InputPort(const std::string &name, const std::map<std::string, std::any> &vars):
            CellBase(1, 1, name) { device_vars = vars; }

        ~InputPort() override = default;

        void update(double wavelength) override { E_out = E_in; }
        void control(std::map<std::string, std::any> ctrl_vars) override {}
    };

    /**
     * Output Port Derived from CellBase
     */
    class OutputPort: public CellBase {
    public:
        OutputPort(const std::string &name, const std::map<std::string, std::any> &vars):
                CellBase(1, 1, name) { device_vars = vars; }

        ~OutputPort() override = default;

        void update(double wavelength) override { E_out = E_in; }
        void control(std::map<std::string, std::any> ctrl_vars) override {}
    };
}

#endif //OMNICLO_TERMINATOR_HPP
