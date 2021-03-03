#!/usr/bin/python3
"""
Module to generate a POTCAR file from n provided files
"""

import argparse as ap
import sys
import numpy as np

class GeneratePOTCAR:
    """
    Class to generate a POTCAR file from given files

    Attributes
    ----------
    files : list
            List of files to be concatenated
    """

    def __init__(self):
        """
        Python constructor
        """
        self.files = []  # files to be concatenated

        self._argparse()  # run argparse

    def _argparse(self):
        """
        Run argparse methods
        """
        parser = ap.ArgumentParser(description="Concatenate several POTCAR files.")
        parser.add_argument('-f', '--files', nargs='+', type=str, required=True, 
                                                help='files to be concatenated')

        args = parser.parse_args()

        self.files = args.files  # update the class

    def _write_file(self):
        """
        Write a concatenated POTCAR file
        """
        
        potcar = open("POTCAR", 'w')
        for item in self.files:
            f = open(item, 'r', errors='ignore')
            potcar.writelines(f.readlines())
            f.close()
        potcar.close()

    def concatenate(self):
        """
        Run the concatenation
        """
        self._write_file()  # read and write the POTCAR files

def main():
    """
    Main function to run the code
    """
    Concatenator = GeneratePOTCAR()  # instantiate the class
    Concatenator.concatenate()       # perform concatenation

if __name__ == "__main__":
    """
    Standard python boilerplate
    """
    main()  # call the main function

