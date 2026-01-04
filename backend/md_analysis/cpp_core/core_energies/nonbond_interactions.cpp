#include "nonbond_interactions.h"
#include "radii_lists.h"
#include <cmath>
#include <algorithm>

const double COULOMB_CONSTANT = 332.06;

double AtomData::distanceTo(const AtomData& other) const {
    double dx = position[0] - other.position[0];
    double dy = position[1] - other.position[1];
    double dz = position[2] - other.position[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<std::shared_ptr<AtomData>> NonbondedInteractions::getNeighborList(
    const std::vector<MoleculeData>& simulationSpace,
    const AtomData& atom,
    double cutoffDistance) {
    
    std::vector<std::shared_ptr<AtomData>> neighbors;
    
    for (const auto& molecule : simulationSpace) {
        for (const auto& neighbor : molecule.atoms) {
            if (atom.id == neighbor.id) continue;
            
            double distance = atom.distanceTo(neighbor);
            if (distance <= cutoffDistance) {
                neighbors.push_back(std::make_shared<AtomData>(neighbor));
            }
        }
    }
    
    return neighbors;
}

std::vector<std::shared_ptr<AtomData>> NonbondedInteractions::getElectrostaticList(
    const std::vector<std::shared_ptr<AtomData>>& neighborList,
    const AtomData& atom,
    double coulombRadius) {
    
    std::vector<std::shared_ptr<AtomData>> electrostaticNeighbors;
    
    for (const auto& neighbor : neighborList) {
        double distance = atom.distanceTo(*neighbor);
        if (distance <= coulombRadius) {
            electrostaticNeighbors.push_back(neighbor);
        }
    }
    
    return electrostaticNeighbors;
}

std::vector<std::shared_ptr<AtomData>> NonbondedInteractions::getVdwList(
    const std::vector<std::shared_ptr<AtomData>>& neighborList,
    const AtomData& atom,
    double vdwRadius) {
    
    std::vector<std::shared_ptr<AtomData>> vdwNeighbors;
    
    for (const auto& neighbor : neighborList) {
        double distance = atom.distanceTo(*neighbor);
        if (distance <= vdwRadius) {
            vdwNeighbors.push_back(neighbor);
        }
    }
    
    return vdwNeighbors;
}

double NonbondedInteractions::calculateLennardJones(
    const AtomData& atom1,
    const AtomData& atom2,
    double sigma,
    double epsilon) {
    
    double r = atom1.distanceTo(atom2);
    
    if (r < 1e-10) return 1e10; 
    
    double r6 = std::pow(r, 6);
    double r12 = r6 * r6;
    double sigma6 = std::pow(sigma, 6);
    double sigma12 = sigma6 * sigma6;
    
    return 4.0 * epsilon * ((sigma12 / r12) - (sigma6 / r6));
}

double NonbondedInteractions::calculateCoulomb(
    const AtomData& atom1,
    const AtomData& atom2,
    double dielectricConstant) {
    
    double r = atom1.distanceTo(atom2);
    
    if (r < 1e-10) return 1e10;
    
    double chargeProduct = atom1.charge * atom2.charge;
    return (COULOMB_CONSTANT * chargeProduct) / (dielectricConstant * r);
}

NonbondedInteractions::NonbondedEnergy NonbondedInteractions::calculateNonbondedEnergy(
    const AtomData& atom1,
    const AtomData& atom2,
    double sigma,
    double epsilon,
    double dielectricConstant) {
    
    NonbondedEnergy energy;
    energy.lennardJones = calculateLennardJones(atom1, atom2, sigma, epsilon);
    energy.coulomb = calculateCoulomb(atom1, atom2, dielectricConstant);
    energy.total = energy.lennardJones + energy.coulomb;
    
    return energy;
}

std::array<double, 3> NonbondedInteractions::calculateLJForce(
    const AtomData& atom1,
    const AtomData& atom2,
    double sigma,
    double epsilon) {
    
    double r = atom1.distanceTo(atom2);
    
    if (r < 1e-10) return {0.0, 0.0, 0.0};
    
    double ux = (atom1.position[0] - atom2.position[0]) / r;
    double uy = (atom1.position[1] - atom2.position[1]) / r;
    double uz = (atom1.position[2] - atom2.position[2]) / r;
    
    double r6 = std::pow(r, 6);
    double r7 = r * r6;
    double r13 = r7 * r6;
    double sigma6 = std::pow(sigma, 6);
    double sigma12 = sigma6 * sigma6;
    
    double dUdr = 4.0 * epsilon * (-12.0 * sigma12 / r13 + 6.0 * sigma6 / r7);
    
    return {
        dUdr * ux,
        dUdr * uy,
        dUdr * uz
    };
}

std::array<double, 3> NonbondedInteractions::calculateCoulombForce(
    const AtomData& atom1,
    const AtomData& atom2,
    double dielectricConstant) {
    
    double r = atom1.distanceTo(atom2);    
    if (r < 1e-10) return {0.0, 0.0, 0.0};

    double ux = (atom1.position[0] - atom2.position[0]) / r;
    double uy = (atom1.position[1] - atom2.position[1]) / r;
    double uz = (atom1.position[2] - atom2.position[2]) / r;
    
    double chargeProduct = atom1.charge * atom2.charge;
    double dUdr = -(COULOMB_CONSTANT * chargeProduct) / (dielectricConstant * r * r);
    
    return {
        dUdr * ux,
        dUdr * uy,
        dUdr * uz
    };
}
