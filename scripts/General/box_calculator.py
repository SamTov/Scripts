"""
A short script to calculate the volume of a cubic or spherical box given information
about the constituents. Similar to the PackMol calculator function.
"""

import numpy as np


class System:
    """
    A class to define the elements of a system and its box size

    Attributes
    -----------
    species_dictionary : dict
        Dictionary of all species in the system, their masses, and their quantity. 
        e.g. {'Na': (22.989, 100), 'Cl': (35.453, 100)}.
    number_of_molecules : int
            Total number of atoms int thr system.
    density : float
            density of the system.
    volume : float
            volume of the system.
    cubic_box_length : float
            length of a cubic box of given volume
    spherical_box_radius : float
            radius of a spherical box of given volume
    """

    def __init__(self):
        """
        Standard constructor
        """
        self.species_dictionary: dict = {}
        self.number_of_molecules: int = 0
        self.density: float = 0.0
        self.volume: float = 0.0

        self.cubic_box_length: float = 0.0
        self.spherical_box_radius: float = 0.0

    def _get_atom_information(self):
        """
        get information from the user about atom types and mass

        Returns
        -------
        updates the self.species_information attribute of the class
        """
        number_of_species: int = int(input("Please enter the number of species in your system: "))

        for i in range(number_of_species):
            symbol: str = input("Species symbol: ")
            mass: str = input("Species mass: ")
            quantity: str = input(f"Number {symbol}s in system: ")
            self.species_dictionary[symbol] = (float(mass), int(quantity))

    def _get_density_information(self):
        """
        get density information of the system

        Returns
        -------
        updates the self.density attribute of the class
        """
        self.density: float = float(input("Density of the system: "))

    def _calculate_volume(self):
        """
        Using class information, calculate the volume of the system

        Returns
        -------
        updates the self.volume attribute of the class
        """
        mass: float = 0.0

        for item in self.species_dictionary:
            mass += self.species_dictionary[item][0] * self.species_dictionary[item][1]

        self.volume = mass / self.density

    def _calculate_cubic_box(self):
        """
        Calculate the cubic box size

        Returns
        -------
        updates the self.cubic_box_length class attribute
        """
        self.cubic_box_length = self.volume ** (1 / 3)

    def _calculate_spherical_box(self):
        """
        Calculate the radius of a spherical system

        Returns
        -------
        updates the self.spherical_box_radius class attribute
        """
        self.spherical_box_radius = (3 * self.volume / (4 * np.pi)) ** (1 / 3)

    def get_box_properties(self):
        """
        Run all arguments and print the box properties.

        Returns
        -------
        Prints output to screen.
        """

        self._get_atom_information()  # collect the species information
        self._get_density_information()  # collect the density information
        self._calculate_volume()  # calculate the volume of the system
        self._calculate_cubic_box()  # calculate length of cubic box
        self._calculate_spherical_box()  # calculate radius of spherical box

        # print results to the screen
        print(f"Cubic box length: {self.cubic_box_length}")
        print(f"Spherical box radius: {self.spherical_box_radius}")


def _main():
    """
    Main function to collect other functions and perform analysis
    """
    system = System()  # instantiate class
    system.get_box_properties()  # run the calculations


if __name__ == "__main__":
    """
    Typical python boilerplate
    """
    _main()
