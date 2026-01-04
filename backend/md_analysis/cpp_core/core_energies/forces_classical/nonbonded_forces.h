#pragma once

#include <vector>
#include <array>
#include <memory>
#include <cmath>
#include <string>

struct AtomData {
    int id;
    std::string element;
    std::array<double, 3> velocity;
    double mass;
    int atomicNumber;
    double charge;
    
    double distanceTo(const AtomData& other) const;
};

struct MoleculeData {
    std::vector<AtomData> atoms;
};

class NonbondedForces {
public: 

    static void computeElectrostaticForces(
        const std::vector<std::vector<AtomData>> neighborLists, 
        const std::vector<std::array<double,3>>& positions,
        const std::vector<double>& charges,
        std::vector<std::array<double,3>>& forces, 
        double radiusCoulombic
    );

    static void computeLJForces(
        const std::vector<std::vector<AtomData>> neighborLists, 
        const std::vector<std::array<double,3>>& positions,
        const std::vector<double>& sigma, 
        const std::vector<double>& epsilon,
        std::vector<std::array<double,3>>& forces, 
        double radiusVDW
    );
};