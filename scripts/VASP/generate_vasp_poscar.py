#!/usr/bin/python3
"""
Script to generate a POSCAR file from an xyz file
"""

import argparse as ap
import sys
import numpy as np

class XYZ2POSCAR:
    """
    Class for converting xyz files into POSCAR files

    Attributes
    ----------
    filepath : str
            Path to they xyz file being converted
    species : dict
            Information about the species in the system.
    """

    def __init__(self):
        """
        Constructor method.
        """

        self.filepath = None     # path to the xyz file
        self.species = {}        # information about the species 
        self.temperature = None  # temperature of the simulation
        self.density = None      # density of the box
        self.box = None          # box vectors of the system

        self._argparse()         # run argparse and update class

    def _argparse(self):
        """
        Run the argparse methods
        """
        
        parser = ap.ArgumentParser(description='Convert xyz files to POSCAR file.')
        parser.add_argument('-f', '--filename', type=str, required=True, 
                                                help='path to xyz file')
        parser.add_argument('-bx', '--box_x', nargs='+', type=float, required=True,
                                                                help='x box vector')
        parser.add_argument('-by', '--box_y', nargs='+', type=float, required=True,
                                                                help='y box vector')
        parser.add_argument('-bz', '--box_z', nargs='+', type=float, required=True,
                                                                help='z box vector')
        parser.add_argument('-t', '--temperature', required=False, type=float,
                                                help='Temperature of  the system')
        parser.add_argument('-p', '--rho', required=False, type=float,
                                        help='Density of the system')

        args = parser.parse_args()  # parse the args

        # update the class
        self.filepath = args.filename
        self.box = [args.box_x, args.box_y, args.box_z]
        if args.temperature:
            self.temperature = args.temperature
        if args.rho:
            self.density = args.rho

    def _read_file(self):
        """
        Read the datafile

        Returns
        -------
        data_array : list
                
        """
        f = open(self.filepath, 'r')
        data_array = []
        
        number_of_atoms = int(f.readline())
        f.readline()  # skip the comment line
        
        counter = 0  # instantiate counter
        for line in f:
            data = line.split()
            if data[0] in self.species:
                self.species[data[0]].append(counter)
            else:
                self.species[data[0]] = []
                self.species[data[0]].append(counter)
            
            data_array.append(data)
            counter += 1
        
        f.close()  # close the file
        test_number = 0
        for item in self.species:
            test_number += len(self.species[item])


        if test_number != number_of_atoms:
            print("Number of atoms in header is not equal to the number of coordinates.")
            sys.exit(1)

        return data_array
    
    def _write_header(self):
        """
        return a header line with the available information

        Returns
        -------
        header : str
                header line to be written in the POSCAR
        """
        
        elements = list(self.species)
        element_string = ""
        for element in elements:
            element_string += f"{element}{len(self.species[element])}"

        return f"POSCAR: {element_string} Temperature: {self.temperature} K Density: {self.density} g/cm^3"
        

    def _write_poscar(self, coordinates):
        """
        Write a poscar file

        Parameters
        ----------
        coordinates : list
                coordinates to ebe written to the POSCAR        
        """
        header = self._write_header()  # get the header line
        f = open('POSCAR', 'w')  # open a POSCAR file object
        f.write(f"{header} \n")
        f.write("1.0 \n")
        for item in self.box:
            f.write(f"{item[0]} {item[1]} {item[2]} \n")

        species = list(self.species)
        species_string = ""
        for item in species:
            species_string += f"{len(self.species[item])} "
        f.write(f"{species_string} \n")
        f.write("Cartesian \n")
        for species in self.species:
            data = np.array(coordinates)[self.species[species]][:, 1:4]
            for line in data:
                f.write(f"{line[0]}  {line[1]}  {line[2]}  \n")
        f.close()

    def convert_file(self):
        """
        Collect methods and perform conversion
        """
        coordinates = self._read_file()  # read the file
        self._write_poscar(coordinates)  # write the poscar

def main():
    """
    Main file to run the code
    """
    Converter = XYZ2POSCAR()  # instantiate the class
    Converter.convert_file()  # run the conversion

if __name__ == "__main__":
    """
    Standard python boilerplate.
    """
    main()  # call the main function

