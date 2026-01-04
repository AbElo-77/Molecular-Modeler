#include "radii_lists.h"
#include <stdexcept>

const std::unordered_map<std::string, double> RadiiLists::VDW_RADII = {
    {"H", 1.20}, {"C", 1.70}, {"N", 1.55}, {"O", 1.52},
    {"F", 1.47}, {"P", 1.80}, {"S", 1.80}, {"Cl", 1.75},
    {"Br", 1.85}, {"I", 1.98}, {"Si", 2.10}, {"Se", 1.90},
    {"As", 1.85}, {"B", 1.92}, {"Li", 1.82}, {"Na", 2.27},
    {"K", 2.75}, {"Ca", 2.31}, {"Zn", 1.39}, {"Cu", 1.40},
    {"Fe", 1.32}, {"Mg", 1.73}
};

const std::unordered_map<std::string, double> RadiiLists::COVALENT_RADII = {
    {"H", 0.31}, {"C", 0.76}, {"N", 0.71}, {"O", 0.66},
    {"F", 0.57}, {"P", 1.07}, {"S", 1.05}, {"Cl", 1.02},
    {"Br", 1.20}, {"I", 1.39}, {"Si", 1.11}, {"Se", 1.19},
    {"As", 1.21}, {"B", 0.84}, {"Li", 1.28}, {"Na", 1.66},
    {"K", 1.96}, {"Ca", 1.71}, {"Zn", 1.22}, {"Cu", 1.32},
    {"Fe", 1.32}, {"Mg", 1.41}
};

const std::unordered_map<std::string, double> RadiiLists::COULOMB_RADII = {
    {"H", 0.50}, {"C", 1.00}, {"N", 1.10}, {"O", 1.20},
    {"F", 1.30}, {"P", 1.40}, {"S", 1.50}, {"Cl", 1.60},
    {"Br", 1.70}, {"I", 1.80}, {"Si", 1.30}, {"Se", 1.60},
    {"As", 1.65}, {"B", 0.95}, {"Li", 0.76}, {"Na", 1.02},
    {"K", 1.38}, {"Ca", 1.14}, {"Zn", 0.87}, {"Cu", 0.96},
    {"Fe", 1.04}, {"Mg", 0.86}
};

const std::unordered_map<std::string, int> RadiiLists::ATOMIC_NUMBERS = {
    {"H", 1}, {"C", 6}, {"N", 7}, {"O", 8},
    {"F", 9}, {"P", 15}, {"S", 16}, {"Cl", 17},
    {"Br", 35}, {"I", 53}, {"Si", 14}, {"Se", 34},
    {"As", 33}, {"B", 5}, {"Li", 3}, {"Na", 11},
    {"K", 19}, {"Ca", 20}, {"Zn", 30}, {"Cu", 29},
    {"Fe", 26}, {"Mg", 12}
};

double RadiiLists::getVdwRadius(const std::string& element) {
    auto it = VDW_RADII.find(element);
    if (it == VDW_RADII.end()) {
        throw std::runtime_error("Element '" + element + "' not found in VDW radii table");
    }
    return it->second;
}

double RadiiLists::getCovalentRadius(const std::string& element) {
    auto it = COVALENT_RADII.find(element);
    if (it == COVALENT_RADII.end()) {
        throw std::runtime_error("Element '" + element + "' not found in covalent radii table");
    }
    return it->second;
}

double RadiiLists::getCoulombRadius(const std::string& element) {
    auto it = COULOMB_RADII.find(element);
    if (it == COULOMB_RADII.end()) {
        throw std::runtime_error("Element '" + element + "' not found in Coulomb radii table");
    }
    return it->second;
}

int RadiiLists::getAtomicNumber(const std::string& element) {
    auto it = ATOMIC_NUMBERS.find(element);
    if (it == ATOMIC_NUMBERS.end()) {
        throw std::runtime_error("Element '" + element + "' not found in atomic numbers table");
    }
    return it->second;
}

bool RadiiLists::hasElement(const std::string& element) {
    return VDW_RADII.find(element) != VDW_RADII.end();
}
