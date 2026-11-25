"""Fuel-agnostic fuel definitions, stoichiometry, and reactant/product gas mixtures.

Any fuel of the form $\\text{C}_x\\text{H}_y\\text{O}_z$ (including pure hydrogen,
$x = z = 0$) is described by its atomic composition and lower heating value.
Stoichiometry follows the global reaction for complete combustion in dry air
(Heywood, 2018, Sec. 3.3; Turns, 2012, Ch. 2):

$$\\text{C}_x\\text{H}_y\\text{O}_z + \\lambda a_s\\,(\\text{O}_2 + \\dots)
  \\rightarrow x\\,\\text{CO}_2 + \\tfrac{y}{2}\\,\\text{H}_2\\text{O}
  + (\\lambda - 1) a_s\\,\\text{O}_2 + \\dots,
  \\qquad a_s = x + \\frac{y}{4} - \\frac{z}{2}$$

Burned-gas compositions assume complete combustion without dissociation, which is
the standard engineering approximation for heat-release analysis below roughly
2200 K (Heywood, 2018, Sec. 4.7). Rich mixtures ($\\lambda < 1$) are supported
only for carbon-free fuels, where the excess fuel ($\\text{H}_2$) is the sole
additional product; rich hydrocarbon products require a water-gas shift
equilibrium and are rejected.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import ArrayLike, NDArray

from combusta.exceptions import PhysicalBoundsError
from combusta.thermo.nasa_polynomials import NASAGasMixture, SPECIES_DB, mix_nasa_species

FloatArray = NDArray[np.float64]

# IUPAC standard atomic weights [kg/mol]
ATOMIC_MASS_C: Final[float] = 12.0107e-3
ATOMIC_MASS_H: Final[float] = 1.00794e-3
ATOMIC_MASS_O: Final[float] = 15.9994e-3

# Dry air composition by mole fraction (same basis as ``DRY_AIR``)
AIR_MOLE_FRACTIONS: Final[Mapping[str, float]] = {
    "N2": 0.78084,
    "O2": 0.20946,
    "AR": 0.00934,
    "CO2": 0.00036,
}


def _air_molar_mass() -> float:
    return sum(y * SPECIES_DB[sp].molar_mass_kg_per_mol for sp, y in AIR_MOLE_FRACTIONS.items())


AIR_MOLAR_MASS_KG_PER_MOL: Final[float] = _air_molar_mass()
MOLES_AIR_PER_MOLE_O2: Final[float] = 1.0 / AIR_MOLE_FRACTIONS["O2"]


@dataclass(frozen=True, slots=True)
class Fuel:
    """Fuel of global composition $\\text{C}_x\\text{H}_y\\text{O}_z$ with its lower heating value.

    Parameters
    ----------
    name : str
        Identifier.
    carbon_atoms, hydrogen_atoms, oxygen_atoms : float
        Atoms per (average) fuel molecule $x, y, z$ [-]. Non-integer values describe
        surrogate or average formulas.
    lhv_j_per_kg : float
        Lower heating value [J/kg].
    vapor_species : str | None, default=None
        Key in ``SPECIES_DB`` for the gaseous fuel. Required to include fuel vapor in a
        premixed unburned charge (e.g. port-injected hydrogen or natural gas).
    reference : str, default=""
        Source of the composition and heating value.
    """

    name: str
    carbon_atoms: float
    hydrogen_atoms: float
    lhv_j_per_kg: float
    oxygen_atoms: float = 0.0
    vapor_species: str | None = None
    reference: str = ""

    def __post_init__(self) -> None:
        if min(self.carbon_atoms, self.hydrogen_atoms, self.oxygen_atoms) < 0.0:
            raise PhysicalBoundsError("Fuel atom counts must be non-negative.")
        if self.carbon_atoms + self.hydrogen_atoms <= 0.0:
            raise PhysicalBoundsError("Fuel must contain carbon and/or hydrogen.")
        if not np.isfinite(self.lhv_j_per_kg) or self.lhv_j_per_kg <= 0.0:
            raise PhysicalBoundsError("Lower heating value must be finite and > 0 J/kg.")
        if self.stoich_o2_moles <= 0.0:
            raise PhysicalBoundsError("Fuel requires a positive stoichiometric oxygen demand.")
        if self.vapor_species is not None and self.vapor_species not in SPECIES_DB:
            raise KeyError(f"vapor_species '{self.vapor_species}' not found in SPECIES_DB.")

    @property
    def molar_mass_kg_per_mol(self) -> float:
        """Fuel molar mass $W_f$ [kg/mol]."""
        return (
            self.carbon_atoms * ATOMIC_MASS_C
            + self.hydrogen_atoms * ATOMIC_MASS_H
            + self.oxygen_atoms * ATOMIC_MASS_O
        )

    @property
    def stoich_o2_moles(self) -> float:
        """Stoichiometric oxygen demand $a_s = x + y/4 - z/2$ [mol O2 / mol fuel]."""
        return self.carbon_atoms + 0.25 * self.hydrogen_atoms - 0.5 * self.oxygen_atoms

    @property
    def stoichiometric_afr(self) -> float:
        """Stoichiometric air-fuel mass ratio $(A/F)_s$ [-]."""
        moles_air = self.stoich_o2_moles * MOLES_AIR_PER_MOLE_O2
        return moles_air * AIR_MOLAR_MASS_KG_PER_MOL / self.molar_mass_kg_per_mol

    @property
    def is_carbon_free(self) -> bool:
        """``True`` for carbon-free fuels such as hydrogen."""
        return self.carbon_atoms == 0.0


HYDROGEN: Final[Fuel] = Fuel(
    name="hydrogen",
    carbon_atoms=0.0,
    hydrogen_atoms=2.0,
    lhv_j_per_kg=119.96e6,
    vapor_species="H2",
    reference="H2; LHV from NIST-JANAF enthalpy of formation of H2O(g)",
)
METHANE: Final[Fuel] = Fuel(
    name="methane",
    carbon_atoms=1.0,
    hydrogen_atoms=4.0,
    lhv_j_per_kg=50.03e6,
    vapor_species="CH4",
    reference="CH4; Heywood (2018) App. D",
)
ISO_OCTANE: Final[Fuel] = Fuel(
    name="iso-octane",
    carbon_atoms=8.0,
    hydrogen_atoms=18.0,
    lhv_j_per_kg=44.3e6,
    reference="C8H18; Heywood (2018) App. D",
)
N_DODECANE: Final[Fuel] = Fuel(
    name="n-dodecane",
    carbon_atoms=12.0,
    hydrogen_atoms=26.0,
    lhv_j_per_kg=44.11e6,
    reference="C12H26; ECN Spray A fuel, LHV from NIST heat of combustion",
)
DIESEL_SURROGATE: Final[Fuel] = Fuel(
    name="diesel (C12H23 surrogate)",
    carbon_atoms=12.0,
    hydrogen_atoms=23.0,
    lhv_j_per_kg=42.8e6,
    reference="Average formula C12H23; typical EN 590 LHV. Override with the fuel certificate.",
)
JET_A1_SURROGATE: Final[Fuel] = Fuel(
    name="Jet A-1 / JP-8 (C12H23 surrogate)",
    carbon_atoms=12.0,
    hydrogen_atoms=23.0,
    lhv_j_per_kg=43.2e6,
    reference="Average formula C12H23 (NASA CEA 'Jet-A'); typical LHV, ASTM D1655 minimum 42.8 MJ/kg.",
)

FUEL_LIBRARY: Final[Mapping[str, Fuel]] = {
    fuel.name: fuel
    for fuel in (HYDROGEN, METHANE, ISO_OCTANE, N_DODECANE, DIESEL_SURROGATE, JET_A1_SURROGATE)
}


def _check_lambda(lambda_air: float) -> float:
    lam = float(lambda_air)
    if not np.isfinite(lam) or lam <= 0.0:
        raise PhysicalBoundsError(f"Air-excess ratio lambda must be finite and > 0, got {lam}.")
    return lam


def reactant_moles(
    fuel: Fuel,
    lambda_air: float,
    *,
    include_fuel_vapor: bool | None = None,
) -> dict[str, float]:
    """Species moles of the unburned charge per mole of fuel [mol].

    Parameters
    ----------
    fuel : Fuel
        Fuel definition.
    lambda_air : float
        Air-excess ratio $\\lambda = (A/F) / (A/F)_s$ [-].
    include_fuel_vapor : bool | None, default=None
        Include gaseous fuel in the charge. ``None`` includes it when the fuel has a
        ``vapor_species`` (premixed gaseous fuel); direct-injected liquid fuels then
        contribute air only during compression.
    """
    lam = _check_lambda(lambda_air)
    moles_air = lam * fuel.stoich_o2_moles * MOLES_AIR_PER_MOLE_O2
    moles = {sp: y * moles_air for sp, y in AIR_MOLE_FRACTIONS.items()}
    include = fuel.vapor_species is not None if include_fuel_vapor is None else include_fuel_vapor
    if include:
        if fuel.vapor_species is None:
            raise PhysicalBoundsError(
                f"Fuel '{fuel.name}' has no vapor_species; cannot include fuel vapor."
            )
        moles[fuel.vapor_species] = moles.get(fuel.vapor_species, 0.0) + 1.0
    return moles


def product_moles(fuel: Fuel, lambda_air: float) -> dict[str, float]:
    """Complete-combustion product moles per mole of fuel [mol] (no dissociation).

    Lean and stoichiometric mixtures ($\\lambda \\ge 1$) are supported for every fuel.
    Rich mixtures ($\\lambda < 1$) are supported for carbon-free fuels, where all
    oxygen forms water and the excess hydrogen remains as $\\text{H}_2$.
    """
    lam = _check_lambda(lambda_air)
    a_s = fuel.stoich_o2_moles
    moles_air = lam * a_s * MOLES_AIR_PER_MOLE_O2
    air = {sp: y * moles_air for sp, y in AIR_MOLE_FRACTIONS.items()}

    products: dict[str, float] = {
        "N2": air["N2"],
        "AR": air["AR"],
        "CO2": air["CO2"] + fuel.carbon_atoms,
    }
    if lam >= 1.0:
        products["H2O"] = 0.5 * fuel.hydrogen_atoms
        products["O2"] = (lam - 1.0) * a_s
        return products
    if not fuel.is_carbon_free:
        raise PhysicalBoundsError(
            "Rich (lambda < 1) products of carbon-containing fuels need a water-gas shift "
            "equilibrium and are not supported."
        )
    o_available = 2.0 * air["O2"] + fuel.oxygen_atoms
    products["H2O"] = o_available
    products["H2"] = 0.5 * fuel.hydrogen_atoms - o_available
    products["O2"] = 0.0
    return products


def element_moles(species_moles: Mapping[str, float], fuel: Fuel | None = None) -> dict[str, float]:
    """Count C, H, O, N, and Ar atoms in a species mixture [mol].

    Fuel molecules without a ``SPECIES_DB`` entry can be counted by passing their
    moles under the key ``"FUEL"`` together with ``fuel``.
    """
    atoms = {
        "N2": {"N": 2.0},
        "O2": {"O": 2.0},
        "AR": {"AR": 1.0},
        "CO2": {"C": 1.0, "O": 2.0},
        "H2O": {"H": 2.0, "O": 1.0},
        "H2": {"H": 2.0},
        "CH4": {"C": 1.0, "H": 4.0},
    }
    totals = {"C": 0.0, "H": 0.0, "O": 0.0, "N": 0.0, "AR": 0.0}
    for species, n in species_moles.items():
        if species == "FUEL":
            if fuel is None:
                raise ValueError("Pass 'fuel' to count atoms of the 'FUEL' pseudo-species.")
            composition = {"C": fuel.carbon_atoms, "H": fuel.hydrogen_atoms, "O": fuel.oxygen_atoms}
        else:
            composition = atoms[species]
        for element, count in composition.items():
            totals[element] += n * count
    return totals


def _nonzero(moles: Mapping[str, float]) -> dict[str, float]:
    return {sp: n for sp, n in moles.items() if n > 0.0}


def make_reactant_mixture(
    fuel: Fuel,
    lambda_air: float,
    *,
    include_fuel_vapor: bool | None = None,
) -> NASAGasMixture:
    """NASA ideal-gas mixture of the unburned charge (air, plus premixed fuel vapor)."""
    moles = reactant_moles(fuel, lambda_air, include_fuel_vapor=include_fuel_vapor)
    return mix_nasa_species(f"reactants_{fuel.name}_lambda_{lambda_air:.3f}", _nonzero(moles))


def make_product_mixture(fuel: Fuel, lambda_air: float) -> NASAGasMixture:
    """NASA ideal-gas mixture of complete-combustion products (no dissociation)."""
    moles = product_moles(fuel, lambda_air)
    return mix_nasa_species(f"products_{fuel.name}_lambda_{lambda_air:.3f}", _nonzero(moles))


@dataclass(frozen=True, slots=True)
class BlendedProperties:
    """Mass-weighted properties of a partially burned single-zone charge.

    Attributes
    ----------
    cp_j_per_kg_k, cv_j_per_kg_k : NDArray[np.float64]
        Specific heats $c_p, c_v$ [J/(kg K)].
    r_specific_j_per_kg_k : NDArray[np.float64]
        Specific gas constant $R$ [J/(kg K)].
    gamma : NDArray[np.float64]
        Specific heat ratio $\\gamma = c_p / c_v$ [-].
    """

    cp_j_per_kg_k: FloatArray
    cv_j_per_kg_k: FloatArray
    r_specific_j_per_kg_k: FloatArray
    gamma: FloatArray


def blend_burned_unburned(
    temperature_k: ArrayLike,
    mass_fraction_burned: ArrayLike,
    unburned: NASAGasMixture,
    burned: NASAGasMixture,
) -> BlendedProperties:
    """Blend unburned and burned gas properties by mass fraction burned $x_b$.

    A single-zone charge is treated as an ideal mixture of unburned and burned gas
    at the common temperature $T$; mass-specific properties mix linearly in mass
    fraction (Gatowski et al., 1984, SAE 841359; Heywood, 2018, Sec. 9.2):

    $$c_p = (1 - x_b) c_{p,u}(T) + x_b c_{p,b}(T), \\qquad
      R = (1 - x_b) R_u + x_b R_b, \\qquad \\gamma = \\frac{c_p}{c_p - R}$$
    """
    t = np.asarray(temperature_k, dtype=np.float64)
    x_b = np.asarray(mass_fraction_burned, dtype=np.float64)
    if not np.all(np.isfinite(x_b)) or np.any((x_b < 0.0) | (x_b > 1.0)):
        raise PhysicalBoundsError("Mass fraction burned x_b must lie within [0, 1].")
    x_b = np.broadcast_to(x_b, t.shape)
    cp = (1.0 - x_b) * unburned.cp(t) + x_b * burned.cp(t)
    r_mix = (1.0 - x_b) * unburned.specific_gas_constant + x_b * burned.specific_gas_constant
    cv = cp - r_mix
    return BlendedProperties(
        cp_j_per_kg_k=np.asarray(cp, dtype=np.float64),
        cv_j_per_kg_k=np.asarray(cv, dtype=np.float64),
        r_specific_j_per_kg_k=np.asarray(r_mix, dtype=np.float64),
        gamma=np.asarray(cp / cv, dtype=np.float64),
    )
