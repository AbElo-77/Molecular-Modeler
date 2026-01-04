#pragma once

#include <vector>
#include <array>
#include <memory>
#include <cmath>

struct BondData {
    int atom1Id;
    int atom2Id;
    int order;  
    double equilibriumLength;
    double forceConstant; 
    bool isRotatable;
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

class BondedInteractions {
public:

    static double calculateBondEnergy(
        double currentLength,
        double equilibriumLength,
        double forceConstant);
    
    static double calculateBondForce(
        double currentLength,
        double equilibriumLength,
        double forceConstant);
    
    static double calculateAngleEnergy(
        double currentAngle,
        double equilibriumAngle,
        double forceConstant);
    
    static double calculateAngleTorque(
        double currentAngle,
        double equilibriumAngle,
        double forceConstant);
    
    static double calculateDihedralEnergy(
        double currentDihedral,
        double periodicity,
        double barrierHeight,
        double phaseOffset);
    
    static double calculateDihedralTorque(
        double currentDihedral,
        double periodicity,
        double barrierHeight,
        double phaseOffset);
    
    static double calculateImproperEnergy(
        double currentDihedral,
        double forceConstant);
    
    static double calculateBondLength(
        const std::array<double, 3>& pos1,
        const std::array<double, 3>& pos2);
    
    static double calculateAngle(
        const std::array<double, 3>& pos1,
        const std::array<double, 3>& pos2,
        const std::array<double, 3>& pos3);
    
    static double calculateDihedral(
        const std::array<double, 3>& pos1,
        const std::array<double, 3>& pos2,
        const std::array<double, 3>& pos3,
        const std::array<double, 3>& pos4);
    
    struct BondedEnergy {
        double bondEnergy;
        double angleEnergy;
        double dihedralEnergy;
        double improperEnergy;
        double total;
    };
    
    static BondedEnergy calculateTotalBondedEnergy(
        const std::vector<BondData>& bonds,
        const std::vector<AngleData>& angles,
        const std::vector<DihedralData>& dihedrals,
        const std::vector<std::array<double, 3>>& positions);
};
