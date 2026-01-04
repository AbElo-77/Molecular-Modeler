#pragma once

struct PairDistance {
    uint8_t graphDistance;
};

struct BondData {
    int atom1Id;
    int atom2Id;
    int order;  
    double equilibriumLength;
    double forceConstant; 
    bool isRotatable;
};

using ExclusionTable = std::unordered_map<uint64_t, PairDistance>;

class Exclusions {
public:

    static ExclusionTable buildExclusionTable(
        int numAtoms,
        const std::vector<BondData>& bonds,
        int maxDepth = 3  
    );
};