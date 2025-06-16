from parapy.gui import display
from parapy.core import *
from wing import WingSurface
from kbeutils import avl
from find_nodes import CodeAster_primitives
import numpy as np


class WingAVLAnalysis(avl.Interface):
    aircraft = Input(in_tree=True)
    case_settings = Input()

    points_number = Input(14)
    rho = Input(1.2)
    Mach = Input(0.7)
    is_mirrored = Input(True)

    @Input
    def check_nodes(self):
        return CodeAster_primitives()

    @Attribute
    def check(self):
        return self.check_nodes.load_primitives

    @Attribute
    def configuration(self):
        return self.aircraft.avl_configuration

    @Part
    def cases(self):
        return avl.Case(quantify=len(self.case_settings),
                        name=self.case_settings[child.index][0],
                        settings=self.case_settings[child.index][1])

    @Attribute
    def print_all_results(self):
        for case_name, result in self.results.items():
            print(f"\n{'=' * 10} Case: {case_name} {'=' * 10}")
            for key, section in result.items():
                print(f"\n>> {key}:")
                if isinstance(section, dict):
                    for subkey, value in section.items():
                        print(f"  {subkey}: {value}")
                else:
                    print(f"  {section}")

    @Attribute
    def print_results_debug(self):
        print(f"Cases in results: {list(self.results.keys())}")
        for case, data in self.results.items():
            print(f"\nCase: {case}")
            print(f"Keys: {data.keys()}")
            print(f"StripForces: {data.get('StripForces')}")

    @Attribute
    def lift_forces(self) -> List[float]:
        # Take the first case
        first_case_name = list(self.results.keys())[0]
        result = self.results[first_case_name]

        # Retrieve StripForces
        strip_forces = result.get('StripForces', {})

        # Check on the right key
        surface_data = strip_forces.get('None', {})
        ccl_list = surface_data.get('c cl', [])

        if not ccl_list:
            print("No 'c cl' in results")
            return []

        # Calculate lift per strip (raw values)
        v = self.Mach * 343
        q = 0.5 * self.rho * v ** 2
        raw = [q * ccl for ccl in ccl_list]

        # ---------------------------------------------------
        # Resample/interpolate to exactly points_number values
        # ---------------------------------------------------
        n_raw = len(raw)
        n_tgt = self.points_number

        if n_raw == 0:
            return []
        if n_raw == n_tgt:
            return raw

        # position grids from 0.0 to 1.0
        x_raw = np.linspace(0.0, 1.0, n_raw)
        x_tgt = np.linspace(0.0, 1.0, n_tgt)

        sampled = np.interp(x_tgt, x_raw, raw)
        return sampled.tolist()[:(self.points_number // 2)] if self.is_mirrored else sampled.tolist()


if __name__ == '__main__':
    # Create wing instance
    wing = WingSurface(label="wing")

    # Define AVL cases (bijv. vaste invalshoek van 5 graden)
    cases = [
        ("alpha_5deg", {'alpha': 5.0}),
    ]

    # Analyse opzetten en tonen in Parapy GUI
    avl_analysis = WingAVLAnalysis(aircraft=wing, case_settings=cases)
    display(avl_analysis)
