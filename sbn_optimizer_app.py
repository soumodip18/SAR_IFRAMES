"""
IRBAS™ Soil Biological Network (SBN) Optimization Engine
Implements multi-objective mixed-integer non-linear optimization for soil consortium design.
"""

from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import minimize
from typing import Dict, List, Tuple, Set


@dataclass(frozen=True)
class SoilProfile:
    ph: float
    ec_ds_m: float
    soc_percent: float
    cec_meq: float
    avail_p_ppm: float
    avail_n_ppm: float
    clay_fraction: float
    moisture_deficit: float  # -psi_m in MPa


@dataclass(frozen=True)
class BiocharProfile:
    feedstock: str
    pyrolysis_temp_c: float
    specific_surface_area_m2_g: float
    mesopore_volume_cm3_g: float
    ph: float
    cec: float


@dataclass
class MicrobeTaxon:
    taxon_id: str
    genus: str
    species: str
    strain: str
    is_pathogen_or_toxigenic: bool
    functional_capacities: Dict[str, float]  # e.g., {'p_sol': 0.8, 'cellulase': 0.9}
    metabolite_production: Dict[str, float]  # e.g., {'gluconate': 0.7, 'organic_c': 0.4}
    metabolite_uptake: Dict[str, float]      # e.g., {'gluconate': 0.0, 'organic_c': 0.8}
    ph_optimum_range: Tuple[float, float]
    ec_tolerance_max: float
    growth_rate_base: float


class SBNOptimizationEngine:
    def __init__(
        self,
        candidate_pool: List[MicrobeTaxon],
        functional_keys: List[str],
        metabolite_keys: List[str],
        redundancy_penalty_coeff: float = 0.35,
        max_consortium_size: int = 5
    ):
        self.raw_pool = candidate_pool
        self.fn_keys = functional_keys
        self.met_keys = metabolite_keys
        self.gamma = redundancy_penalty_coeff
        self.k_max = max_consortium_size

    def _prune_infeasible_taxa(
        self, 
        soil: SoilProfile, 
        biochar: BiocharProfile
    ) -> List[MicrobeTaxon]:
        """Apply strict biosafety, pH, and osmotic stress pruning."""
        feasible = []
        effective_ph = 0.8 * soil.ph + 0.2 * biochar.ph  # Buffered boundary layer estimate
        
        for m in self.raw_pool:
            # Biosafety absolute override
            if m.is_pathogen_or_toxigenic:
                continue
            # Physiological tolerance gate
            if not (m.ph_optimum_range[0] <= effective_ph <= m.ph_optimum_range[1]):
                continue
            if soil.ec_ds_m > m.ec_tolerance_max:
                continue
            feasible.append(m)
        return feasible

    def _compute_soil_deficits(self, soil: SoilProfile) -> np.ndarray:
        """Compute severity weights for target biological functions."""
        deficits = []
        for fn in self.fn_keys:
            if fn == "p_sol":
                # High deficit if available P is low
                w = np.clip((25.0 - soil.avail_p_ppm) / 25.0, 0.0, 1.0)
            elif fn == "cellulase":
                # High deficit if SOC is low (requires carbon cycling acceleration)
                w = np.clip((2.5 - soil.soc_percent) / 2.5, 0.0, 1.0)
            elif fn == "osmoprotection":
                w = np.clip(soil.ec_ds_m / 8.0, 0.0, 1.0)
            else:
                w = 0.5
            deficits.append(w)
        return np.array(deficits)

    def _build_interaction_matrices(
        self, 
        taxa: List[MicrobeTaxon], 
        biochar: BiocharProfile
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        n = len(taxa)
        cross_feed = np.zeros((n, n))
        antagonism = np.zeros((n, n))
        niche_overlap = np.zeros((n, n))

        # Biochar mesopore spatial diffusion modifier (0.1 to 1.0)
        eta = np.clip(biochar.mesopore_volume_cm3_g / 0.5, 0.1, 1.0)

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                
                # Cross-feeding: donor i -> recipient j
                flux = 0.0
                for met in self.met_keys:
                    prod = taxa[i].metabolite_production.get(met, 0.0)
                    uptk = taxa[j].metabolite_uptake.get(met, 0.0)
                    flux += min(prod, uptk)
                cross_feed[i, j] = eta * flux

                # Antagonism: fungal-bacterial interference or competitive exclusion
                if taxa[i].genus == "Trichoderma" and taxa[j].genus == "Bacillus":
                    antagonism[i, j] = 0.4
                elif taxa[i].genus == "Pseudomonas" and taxa[j].genus == "Aspergillus":
                    antagonism[i, j] = 0.2

                # Niche overlap (Jaccard similarity on functional capacities)
                f_i = np.array([taxa[i].functional_capacities.get(k, 0.0) for k in self.fn_keys])
                f_j = np.array([taxa[j].functional_capacities.get(k, 0.0) for k in self.fn_keys])
                denom = np.sum(np.maximum(f_i, f_j))
                if denom > 1e-6:
                    niche_overlap[i, j] = np.sum(np.minimum(f_i, f_j)) / denom

        return cross_feed, antagonism, niche_overlap

    def _evaluate_stability(
        self, 
        x: np.ndarray, 
        taxa: List[MicrobeTaxon], 
        cross_feed: np.ndarray, 
        antagonism: np.ndarray,
        stress_sigma: np.ndarray
    ) -> float:
        """Evaluate asymptotic stability via spectral radius of the effective Jacobian."""
        active_idx = np.where(x > 1e-3)[0]
        if len(active_idx) == 0:
            return -100.0
        
        k = len(active_idx)
        J = np.zeros((k, k))
        
        for idx_i, i in enumerate(active_idx):
            for idx_j, j in enumerate(active_idx):
                if idx_i == idx_j:
                    # Self-limiting carrying capacity under perturbation
                    stress_decay = np.sum(stress_sigma) * 0.1
                    J[idx_i, idx_j] = -1.0 * (taxa[i].growth_rate_base + stress_decay) * x[i]
                else:
                    interaction = cross_feed[i, j] - antagonism[i, j]
                    J[idx_i, idx_j] = interaction * x[i]
                    
        eigenvalues = np.linalg.eigvals(J)
        max_real = np.max(np.real(eigenvalues))
        
        # Stability reward: negative max real eigenvalue indicates return to equilibrium
        if max_real < -1e-4:
            return 1.0 / np.abs(max_real)
        else:
            return -50.0  # Penalty for unstable/explosive configurations

    def optimize(
        self, 
        soil: SoilProfile, 
        biochar: BiocharProfile, 
        stress: np.ndarray
    ) -> Dict:
        feasible_taxa = self._prune_infeasible_taxa(soil, biochar)
        n = len(feasible_taxa)
        
        if n == 0:
            return {"status": "FAILED", "reason": "No candidate survived pruning filters."}

        deficits = self._compute_soil_deficits(soil)
        C_mat, A_mat, O_mat = self._build_interaction_matrices(feasible_taxa, biochar)
        
        # Assemble functional matrix
        M_F = np.zeros((n, len(self.fn_keys)))
        for i, t in enumerate(feasible_taxa):
            for k_idx, k in enumerate(self.fn_keys):
                M_F[i, k_idx] = t.functional_capacities.get(k, 0.0)

        # Objective function for continuous abundance allocation (fixed binary subset or relaxation)
        def objective(x: np.ndarray) -> float:
            z = (x > 1e-3).astype(float)
            
            # 1. Functional Coverage
            covered = 1.0 - np.exp(-np.dot(x, M_F))
            F_score = np.sum(deficits * covered) * 10.0
            
            # 2. Cross-feeding Flux
            C_score = np.sum(np.outer(x, x) * C_mat) * 5.0
            
            # 3. Antagonism & Redundancy Penalty
            A_score = np.sum(np.outer(x, x) * A_mat) * 8.0 + self.gamma * np.sum(np.outer(np.sqrt(x), np.sqrt(x)) * O_mat)
            
            # 4. Mycorrhizal Bridging Module
            has_amf = any(t.genus == "Glomus" for t, xi in zip(feasible_taxa, x) if xi > 1e-3)
            has_psb = any(t.functional_capacities.get("p_sol", 0) > 0.5 for t, xi in zip(feasible_taxa, x) if xi > 1e-3)
            Phi_score = 4.0 if (has_amf and has_psb) else 0.0
            
            # 5. Dynamic Stability
            R_score = self._evaluate_stability(x, feasible_taxa, C_mat, A_mat, stress)
            
            # Sparsity penalty (soft constraint on consortium size)
            size_penalty = 15.0 * max(0, np.sum(z) - self.k_max)
            
            total_sbn = F_score + C_score + Phi_score + R_score - A_score - size_penalty
            return -total_sbn  # Minimize negative objective

        # Optimize continuous titers x in [0, 1]
        x0 = np.ones(n) * 0.2
        bounds = [(0.0, 1.0) for _ in range(n)]
        res = minimize(objective, x0, bounds=bounds, method="SLSQP")

        best_x = res.x
        selected_consortium = []
        for idx, abundance in enumerate(best_x):
            if abundance > 0.05:
                selected_consortium.append({
                    "taxon": f"{feasible_taxa[idx].genus} {feasible_taxa[idx].species} ({feasible_taxa[idx].strain})",
                    "titer_fraction": round(float(abundance), 3),
                    "primary_role": max(feasible_taxa[idx].functional_capacities, key=feasible_taxa[idx].functional_capacities.get)
                })

        return {
            "status": "SUCCESS",
            "sbn_score": round(float(-res.fun), 4),
            "consortium_size": len(selected_consortium),
            "recommended_consortium": selected_consortium,
            "boundary_conditions": {
                "biochar_feedstock": biochar.feedstock,
                "pyrolysis_temp": biochar.pyrolysis_temp_c,
                "effective_ph": round(0.8 * soil.ph + 0.2 * biochar.ph, 2)
            }
        }


# =====================================================================
# Verification Execution
# =====================================================================
if __name__ == "__main__":
    # Test Data: Degraded Saline Soil with Phosphorus Lockup
    soil_sample = SoilProfile(
        ph=8.2,
        ec_ds_m=4.5,
        soc_percent=0.6,
        cec_meq=14.0,
        avail_p_ppm=6.2,
        avail_n_ppm=18.0,
        clay_fraction=0.35,
        moisture_deficit=1.2
    )

    biochar_sample = BiocharProfile(
        feedstock="Cotton_Stalk_Slow_Pyrolysis",
        pyrolysis_temp_c=550.0,
        specific_surface_area_m2_g=380.0,
        mesopore_volume_cm3_g=0.32,
        ph=7.8,
        cec=24.0
    )

    # Candidate Microbial Library
    library = [
        MicrobeTaxon("M01", "Aspergillus", "fumigatus", "AF-01", True, {"cellulase": 0.9}, {"organic_c": 0.8}, {}, (5.0, 8.0), 5.0, 0.4),
        MicrobeTaxon("M02", "Aspergillus", "niger", "AN-MBF7", False, {"p_sol": 0.95, "cellulase": 0.7}, {"gluconate": 0.9, "organic_c": 0.6}, {}, (5.5, 8.8), 6.0, 0.35),
        MicrobeTaxon("M03", "Bacillus", "subtilis", "BS-MBF2", False, {"osmoprotection": 0.8, "eps": 0.85}, {"phytohormone": 0.6}, {"gluconate": 0.7, "organic_c": 0.5}, (6.0, 9.0), 7.0, 0.5),
        MicrobeTaxon("M04", "Penicillium", "bilaiae", "PB-04", False, {"p_sol": 0.85}, {"citrate": 0.7}, {}, (5.0, 7.5), 3.0, 0.3),
        MicrobeTaxon("M05", "Pseudomonas", "putida", "PP-MBF1", False, {"siderophore": 0.9, "osmoprotection": 0.6}, {}, {"organic_c": 0.7}, (6.0, 8.5), 4.0, 0.45),
        MicrobeTaxon("M06", "Glomus", "intraradices", "AMF-01", False, {"mycorrhizal_bridge": 1.0}, {}, {}, (6.0, 8.5), 5.0, 0.2)
    ]

    engine = SBNOptimizationEngine(
        candidate_pool=library,
        functional_keys=["p_sol", "cellulase", "osmoprotection"],
        metabolite_keys=["gluconate", "organic_c", "citrate"]
    )

    stress_perturbation = np.array([1.5, 0.8, 0.5])  # Heat, drought, osmotic increments
    result = engine.optimize(soil_sample, biochar_sample, stress_perturbation)
    
    print("Optimization Execution Output:")
    import pprint
    pprint.pprint(result)
