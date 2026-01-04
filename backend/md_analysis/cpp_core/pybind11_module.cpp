#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "radii_lists.h"
#include "nonbond_interactions.h"
#include "bonded_interactions.h"

namespace py = pybind11;

PYBIND11_MODULE(molecular_interactions, m) {
    m.doc() = "C++ molecular interaction calculations with pybind11";
    
    py::class_<RadiiLists>(m, "RadiiLists", "Van der Waals and interaction radii data")
        .def_static("get_vdw_radius", &RadiiLists::getVdwRadius, 
                   "Get van der Waals radius for an element")
        .def_static("get_covalent_radius", &RadiiLists::getCovalentRadius,
                   "Get covalent radius for an element")
        .def_static("get_coulomb_radius", &RadiiLists::getCoulombRadius,
                   "Get Coulomb interaction radius for an element")
        .def_static("get_atomic_number", &RadiiLists::getAtomicNumber,
                   "Get atomic number for an element")
        .def_static("has_element", &RadiiLists::hasElement,
                   "Check if element exists in radii tables")
        .def_static("get_vdw_radii_dict", []() {
            py::dict d;
            for (const auto& pair : RadiiLists::VDW_RADII) {
                d[pair.first.c_str()] = pair.second;
            }
            return d;
        }, "Get all VDW radii as dictionary")
        .def_static("get_covalent_radii_dict", []() {
            py::dict d;
            for (const auto& pair : RadiiLists::COVALENT_RADII) {
                d[pair.first.c_str()] = pair.second;
            }
            return d;
        }, "Get all covalent radii as dictionary")
        .def_static("get_coulomb_radii_dict", []() {
            py::dict d;
            for (const auto& pair : RadiiLists::COULOMB_RADII) {
                d[pair.first.c_str()] = pair.second;
            }
            return d;
        }, "Get all Coulomb radii as dictionary");
    
    py::class_<AtomData>(m, "AtomData", "Atom data structure")
        .def(py::init<>())
        .def_readwrite("id", &AtomData::id)
        .def_readwrite("element", &AtomData::element)
        .def_readwrite("position", &AtomData::position)
        .def_readwrite("velocity", &AtomData::velocity)
        .def_readwrite("mass", &AtomData::mass)
        .def_readwrite("atomic_number", &AtomData::atomicNumber)
        .def_readwrite("charge", &AtomData::charge)
        .def("distance_to", &AtomData::distanceTo, "Calculate distance to another atom");
    
    py::class_<MoleculeData>(m, "MoleculeData", "Molecule data structure")
        .def(py::init<>())
        .def_readwrite("atoms", &MoleculeData::atoms);
    
    py::class_<NonbondedInteractions>(m, "NonbondedInteractions", 
                                      "Nonbonded interaction calculations")
        .def_static("get_neighbor_list", &NonbondedInteractions::getNeighborList,
                   py::arg("simulation_space"),
                   py::arg("atom"),
                   py::arg("cutoff_distance") = 0.55,
                   "Get neighbor list for an atom")
        .def_static("get_electrostatic_list", &NonbondedInteractions::getElectrostaticList,
                   "Get electrostatic neighbors")
        .def_static("get_vdw_list", &NonbondedInteractions::getVdwList,
                   "Get van der Waals neighbors")
        .def_static("calculate_lennard_jones", &NonbondedInteractions::calculateLennardJones,
                   py::arg("atom1"),
                   py::arg("atom2"),
                   py::arg("sigma"),
                   py::arg("epsilon"),
                   "Calculate Lennard-Jones potential")
        .def_static("calculate_coulomb", &NonbondedInteractions::calculateCoulomb,
                   py::arg("atom1"),
                   py::arg("atom2"),
                   py::arg("dielectric_constant") = 1.0,
                   "Calculate Coulomb potential")
        .def_static("calculate_lj_force", &NonbondedInteractions::calculateLJForce,
                   "Calculate Lennard-Jones force")
        .def_static("calculate_coulomb_force", &NonbondedInteractions::calculateCoulombForce,
                   py::arg("atom1"),
                   py::arg("atom2"),
                   py::arg("dielectric_constant") = 1.0,
                   "Calculate Coulomb force");
    
    py::class_<NonbondedInteractions::NonbondedEnergy>(m, "NonbondedEnergy", 
                                                       "Nonbonded energy components")
        .def(py::init<>())
        .def_readwrite("lennard_jones", &NonbondedInteractions::NonbondedEnergy::lennardJones)
        .def_readwrite("coulomb", &NonbondedInteractions::NonbondedEnergy::coulomb)
        .def_readwrite("total", &NonbondedInteractions::NonbondedEnergy::total);
    
    m.def("calculate_nonbonded_energy", &NonbondedInteractions::calculateNonbondedEnergy,
         py::arg("atom1"),
         py::arg("atom2"),
         py::arg("sigma"),
         py::arg("epsilon"),
         py::arg("dielectric_constant") = 1.0,
         "Calculate combined nonbonded energy");
    
    py::class_<BondData>(m, "BondData", "Bond data structure")
        .def(py::init<>())
        .def_readwrite("atom1_id", &BondData::atom1Id)
        .def_readwrite("atom2_id", &BondData::atom2Id)
        .def_readwrite("order", &BondData::order)
        .def_readwrite("equilibrium_length", &BondData::equilibriumLength)
        .def_readwrite("force_constant", &BondData::forceConstant)
        .def_readwrite("is_rotatable", &BondData::isRotatable);
    
    py::class_<AngleData>(m, "AngleData", "Angle data structure")
        .def(py::init<>())
        .def_readwrite("atom1_id", &AngleData::atom1Id)
        .def_readwrite("atom2_id", &AngleData::atom2Id)
        .def_readwrite("atom3_id", &AngleData::atom3Id)
        .def_readwrite("equilibrium_angle", &AngleData::equilibriumAngle)
        .def_readwrite("force_constant", &AngleData::forceConstant);
    
    py::class_<DihedralData>(m, "DihedralData", "Dihedral/Torsion data structure")
        .def(py::init<>())
        .def_readwrite("atom1_id", &DihedralData::atom1Id)
        .def_readwrite("atom2_id", &DihedralData::atom2Id)
        .def_readwrite("atom3_id", &DihedralData::atom3Id)
        .def_readwrite("atom4_id", &DihedralData::atom4Id)
        .def_readwrite("periodicity", &DihedralData::periodicity)
        .def_readwrite("barrier_height", &DihedralData::barrierHeight)
        .def_readwrite("phase_offset", &DihedralData::phaseOffset);
    
    py::class_<BondedInteractions>(m, "BondedInteractions",
                                   "Bonded interaction calculations")
        .def_static("calculate_bond_energy", &BondedInteractions::calculateBondEnergy,
                   "Calculate bond energy")
        .def_static("calculate_bond_force", &BondedInteractions::calculateBondForce,
                   "Calculate bond force")
        .def_static("calculate_angle_energy", &BondedInteractions::calculateAngleEnergy,
                   "Calculate angle energy")
        .def_static("calculate_angle_torque", &BondedInteractions::calculateAngleTorque,
                   "Calculate angle torque")
        .def_static("calculate_dihedral_energy", &BondedInteractions::calculateDihedralEnergy,
                   "Calculate dihedral energy")
        .def_static("calculate_dihedral_torque", &BondedInteractions::calculateDihedralTorque,
                   "Calculate dihedral torque")
        .def_static("calculate_improper_energy", &BondedInteractions::calculateImproperEnergy,
                   "Calculate improper dihedral energy")
        .def_static("calculate_bond_length", 
                   (double (*)(const std::array<double, 3>&, const std::array<double, 3>&))
                   &BondedInteractions::calculateBondLength,
                   "Calculate bond length between two atoms")
        .def_static("calculate_angle", 
                   (double (*)(const std::array<double, 3>&, 
                               const std::array<double, 3>&, 
                               const std::array<double, 3>&))
                   &BondedInteractions::calculateAngle,
                   "Calculate angle between three atoms")
        .def_static("calculate_dihedral",
                   (double (*)(const std::array<double, 3>&, 
                               const std::array<double, 3>&,
                               const std::array<double, 3>&,
                               const std::array<double, 3>&))
                   &BondedInteractions::calculateDihedral,
                   "Calculate dihedral angle between four atoms");
    
    py::class_<BondedInteractions::BondedEnergy>(m, "BondedEnergy",
                                                 "Bonded energy components")
        .def(py::init<>())
        .def_readwrite("bond_energy", &BondedInteractions::BondedEnergy::bondEnergy)
        .def_readwrite("angle_energy", &BondedInteractions::BondedEnergy::angleEnergy)
        .def_readwrite("dihedral_energy", &BondedInteractions::BondedEnergy::dihedralEnergy)
        .def_readwrite("improper_energy", &BondedInteractions::BondedEnergy::improperEnergy)
        .def_readwrite("total", &BondedInteractions::BondedEnergy::total);
    
    m.def("calculate_total_bonded_energy", &BondedInteractions::calculateTotalBondedEnergy,
         py::arg("bonds"),
         py::arg("angles"),
         py::arg("dihedrals"),
         py::arg("positions"),
         "Calculate total bonded energy");
}
