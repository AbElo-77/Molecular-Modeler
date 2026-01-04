#pragma once

#include <unordered_map>
#include <string>
#include <vector>

class RadiiLists {
public:

    static const std::unordered_map<std::string, double> VDW_RADII;

    static const std::unordered_map<std::string, double> COVALENT_RADII;
    
    static const std::unordered_map<std::string, double> COULOMB_RADII;
    
    static const std::unordered_map<std::string, int> ATOMIC_NUMBERS;
    
    static double getVdwRadius(const std::string& element);
    
    static double getCovalentRadius(const std::string& element);
    
    static double getCoulombRadius(const std::string& element);
    
    static int getAtomicNumber(const std::string& element);
    
    static bool hasElement(const std::string& element);
};
