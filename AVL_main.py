from parapy.gui import display
from parapy.core import *
from wing import WingSurface
from kbeutils import avl


class WingAVLAnalysis(avl.Interface):
    aircraft = Input(in_tree=True)
    case_settings = Input()

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
    def lift_forces(self):
        # Pak eerste case
        first_case_name = list(self.results.keys())[0]
        result = self.results[first_case_name]

        # Haal StripForces op
        strip_forces = result.get('StripForces', {})

        # Check op juiste sleutel
        surface_data = strip_forces.get('None', {})
        ccl_list = surface_data.get('c cl', [])

        if not ccl_list:
            print("No 'c cl' in results")
            return []

        # Bereken lift per strip
        rho = 1200
        v = 0.4 * 343
        q = 0.5 * rho * v ** 2

        return [q * ccl for ccl in ccl_list]


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
