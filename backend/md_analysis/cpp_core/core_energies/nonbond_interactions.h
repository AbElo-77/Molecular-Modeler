#pragma once

#include <vector>
#include <array>
#include <memory>
#include <cmath>

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

struct MoleculeData {
    std::vector<AtomData> atoms;
};

class NonbondedInteractions {
public:

    static std::vector<std::shared_ptr<AtomData>> getNeighborList(
        const std::vector<MoleculeData>& simulationSpace,
        const AtomData& atom,
        double cutoffDistance = 0.55);
    
    static std::vector<std::shared_ptr<AtomData>> getElectrostaticList(
        const std::vector<std::shared_ptr<AtomData>>& neighborList,
        const AtomData& atom,
        double coulombRadius);
    
    static std::vector<std::shared_ptr<AtomData>> getVdwList(
        const std::vector<std::shared_ptr<AtomData>>& neighborList,
        const AtomData& atom,
        double vdwRadius);
    
    static double calculateLennardJones(
        const AtomData& atom1,
        const AtomData& atom2,
        double sigma,
        double epsilon
    );
    
    static double calculateCoulomb(
        const AtomData& atom1,
        const AtomData& atom2,
        double dielectricConstant = 1.0);
    
    struct NonbondedEnergy {
        double lennardJones;
        double coulomb;
        double total;
    };
    
    static NonbondedEnergy calculateNonbondedEnergy(
        const AtomData& atom1,
        const AtomData& atom2,
        double sigma,
        double epsilon,
        double dielectricConstant = 1.0);
    
    static std::array<double, 3> calculateLJForce(
        const AtomData& atom1,
        const AtomData& atom2,
        double sigma,
        double epsilon);
    
    static std::array<double, 3> calculateCoulombForce(
        const AtomData& atom1,
        const AtomData& atom2,
        double dielectricConstant = 1.0);
};
