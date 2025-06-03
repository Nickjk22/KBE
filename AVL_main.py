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
