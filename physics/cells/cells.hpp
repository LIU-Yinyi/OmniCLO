/**
 * @file cells.hpp
 * @brief Integral cells in silicon photonic
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_CELLS_HPP
#define OMNICLO_CELLS_HPP

#include "waveguide.hpp"
#include "coupler.hpp"
#include "microring.hpp"

namespace cells {

    /**
     * Acquire cell type name of which class should be derived from CellBase
     * @tparam T: type of cell
     * @return cell name in string
     */
    template <typename T>
    std::string cell_type_name() {
        if(std::is_same<T, Waveguide>::value) { return "Waveguide"; }
        else if(std::is_same<T, Coupler>::value) { return "Coupler"; }
        else if(std::is_same<T, MicroRingSerial>::value) { return "MicroRingSerial"; }
        else { return ""; }
    }

}

#endif //OMNICLO_CELLS_HPP
