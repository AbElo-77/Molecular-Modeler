#pragma once

#include <vector>
#include <array>

struct AtomData {
    int id;
    std::string element;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    double mass;
    int atomicNumber;
    double charge;
    
    double distanceTo(const AtomData& other) const;
};

struct BondData {
    int atom1Id;
    int atom2Id;
    double equilibriumLength;
    double forceConstant;
};

struct AngleData {
    int atom1Id;
    int atom2Id;
    int atom3Id;
    double equilibriumAngle;
    double forceConstant;
};

struct DihedralData {
    int atom1Id;
    int atom2Id;
    int atom3Id;
    int atom4Id;
    double periodicity;
    double barrierHeight;
    double phaseOffset;
};

class BondedForces {
public:

    static void computeBondForces(
        const std::vector<BondData>& bonds,
        const std::vector<std::array<double,3>>& positions,
        std::vector<std::array<double,3>>& forces
    );

    static void computeAngleForces(
        const std::vector<AngleData>& angles,
        const std::vector<std::array<double,3>>& positions,
        std::vector<std::array<double,3>>& forces
    );

    static void computeDihedralForces(
        const std::vector<DihedralData>& dihedrals,
        const std::vector<std::array<double,3>>& positions,
        std::vector<std::array<double,3>>& forces
    );
};
