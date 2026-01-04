#include "bonded_interactions.h"
#include <cmath>
#include <stdexcept>

const double PI = 3.14159265358979323846;

double BondedInteractions::calculateBondEnergy(
    double currentLength,
    double equilibriumLength,
    double forceConstant) {
    
    double delta = currentLength - equilibriumLength;
    return 0.5 * forceConstant * delta * delta;
}

double BondedInteractions::calculateBondForce(
    double currentLength,
    double equilibriumLength,
    double forceConstant) {
    
    double delta = currentLength - equilibriumLength;
    return forceConstant * delta;
}

double BondedInteractions::calculateAngleEnergy(
    double currentAngle,
    double equilibriumAngle,
    double forceConstant) {
    
    double delta = currentAngle - equilibriumAngle;
    return 0.5 * forceConstant * delta * delta;
}

double BondedInteractions::calculateAngleTorque(
    double currentAngle,
    double equilibriumAngle,
    double forceConstant) {
    
    double delta = currentAngle - equilibriumAngle;
    return forceConstant * delta;
}

double BondedInteractions::calculateDihedralEnergy(
    double currentDihedral,
    double periodicity,
    double barrierHeight,
    double phaseOffset) {
    
    return barrierHeight * (1.0 + std::cos(periodicity * currentDihedral - phaseOffset));
}

double BondedInteractions::calculateDihedralTorque(
    double currentDihedral,
    double periodicity,
    double barrierHeight,
    double phaseOffset) {
    
    return -periodicity * barrierHeight * std::sin(periodicity * currentDihedral - phaseOffset);
}

double BondedInteractions::calculateImproperEnergy(
    double currentDihedral,
    double forceConstant) {
    
    return 0.5 * forceConstant * currentDihedral * currentDihedral;
}

double BondedInteractions::calculateBondLength(
    const std::array<double, 3>& pos1,
    const std::array<double, 3>& pos2) {
    
    double dx = pos1[0] - pos2[0];
    double dy = pos1[1] - pos2[1];
    double dz = pos1[2] - pos2[2];
    
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double BondedInteractions::calculateAngle(
    const std::array<double, 3>& pos1,
    const std::array<double, 3>& pos2,
    const std::array<double, 3>& pos3) {
    
    double v1x = pos1[0] - pos2[0];
    double v1y = pos1[1] - pos2[1];
    double v1z = pos1[2] - pos2[2];
    
    double v2x = pos3[0] - pos2[0];
    double v2y = pos3[1] - pos2[1];
    double v2z = pos3[2] - pos2[2];
    
    double v1_len = std::sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
    double v2_len = std::sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
    
    if (v1_len < 1e-10 || v2_len < 1e-10) {
        throw std::runtime_error("Collinear atoms in angle calculation");
    }
    
    v1x /= v1_len;
    v1y /= v1_len;
    v1z /= v1_len;
    
    v2x /= v2_len;
    v2y /= v2_len;
    v2z /= v2_len;

    double dotProduct = v1x * v2x + v1y * v2y + v1z * v2z;   
    dotProduct = (dotProduct < -1.0) ? -1.0 : (dotProduct > 1.0) ? 1.0 : dotProduct;
    
    return std::acos(dotProduct);
}

double BondedInteractions::calculateDihedral(
    const std::array<double, 3>& pos1,
    const std::array<double, 3>& pos2,
    const std::array<double, 3>& pos3,
    const std::array<double, 3>& pos4) {
    
    double v1x = pos1[0] - pos2[0];
    double v1y = pos1[1] - pos2[1];
    double v1z = pos1[2] - pos2[2];
    
    double v2x = pos3[0] - pos2[0];
    double v2y = pos3[1] - pos2[1];
    double v2z = pos3[2] - pos2[2];
    
    double v3x = pos4[0] - pos3[0];
    double v3y = pos4[1] - pos3[1];
    double v3z = pos4[2] - pos3[2];
    
    double n1x = v1y * v2z - v1z * v2y;
    double n1y = v1z * v2x - v1x * v2z;
    double n1z = v1x * v2y - v1y * v2x;
    
    double n2x = v2y * v3z - v2z * v3y;
    double n2y = v2z * v3x - v2x * v3z;
    double n2z = v2x * v3y - v2y * v3x;
    
    double n1_len = std::sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    double n2_len = std::sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    
    if (n1_len < 1e-10 || n2_len < 1e-10) {
        return 0.0;
    }
    
    n1x /= n1_len;
    n1y /= n1_len;
    n1z /= n1_len;
    
    n2x /= n2_len;
    n2y /= n2_len;
    n2z /= n2_len;
    
    double dotProduct = n1x * n2x + n1y * n2y + n1z * n2z;
    dotProduct = (dotProduct < -1.0) ? -1.0 : (dotProduct > 1.0) ? 1.0 : dotProduct;
    
    double angle = std::acos(dotProduct);
    
    double signx = n1y * n2z - n1z * n2y;
    double signy = n1z * n2x - n1x * n2z;
    double signz = n1x * n2y - n1y * n2x;
    
    double sign = v2x * signx + v2y * signy + v2z * signz;
    if (sign < 0) {
        angle = -angle;
    }
    
    return angle;
}

BondedInteractions::BondedEnergy BondedInteractions::calculateTotalBondedEnergy(
    const std::vector<BondData>& bonds,
    const std::vector<AngleData>& angles,
    const std::vector<DihedralData>& dihedrals,
    const std::vector<std::array<double, 3>>& positions) {
    
    BondedEnergy totalEnergy = {0.0, 0.0, 0.0, 0.0, 0.0};
    
    for (const auto& bond : bonds) {
        double length = calculateBondLength(
            positions[bond.atom1Id],
            positions[bond.atom2Id]);
        totalEnergy.bondEnergy += calculateBondEnergy(
            length,
            bond.equilibriumLength,
            bond.forceConstant);
    }
    
    for (const auto& angle : angles) {
        double angleRad = calculateAngle(
            positions[angle.atom1Id],
            positions[angle.atom2Id],
            positions[angle.atom3Id]);
        totalEnergy.angleEnergy += calculateAngleEnergy(
            angleRad,
            angle.equilibriumAngle,
            angle.forceConstant);
    }
    
    for (const auto& dihedral : dihedrals) {
        double dihedralAngle = calculateDihedral(
            positions[dihedral.atom1Id],
            positions[dihedral.atom2Id],
            positions[dihedral.atom3Id],
            positions[dihedral.atom4Id]);
        totalEnergy.dihedralEnergy += calculateDihedralEnergy(
            dihedralAngle,
            dihedral.periodicity,
            dihedral.barrierHeight,
            dihedral.phaseOffset);
    }
    
    totalEnergy.total = totalEnergy.bondEnergy + 
                       totalEnergy.angleEnergy + 
                       totalEnergy.dihedralEnergy + 
                       totalEnergy.improperEnergy;
    
    return totalEnergy;
}
