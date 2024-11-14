#! /usr/bin/env python3
from os import environ
from sys import argv
from os.path import join
import re

ATOMIC_SYMBOLS = ['neut', 'h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'nh', 'fl', 'mc', 'lv', 'ts', 'og']

def z_to_symbol(z):
    return ATOMIC_SYMBOLS[z]

def symbol_to_z(symbol):
    return ATOMIC_SYMBOLS.index(symbol)

NETS_DIR = join(environ['MESA_DIR'], 'data', 'net_data', 'nets')
def isos_from_net(net_name):
    """Return a list of isotopes from a net file
    
    Parameters
    ----------
    net_name : str
        The name of the net file to read, with or without the .net extension
        
    Returns
    -------
    list of str
        A list of isotopes in the net file, sorted by atomic number and mass number
    """
    # remove the .net extension; will try with and without it later
    net_name = net_name.replace('.net', '')
    # normal nets live in a file with the form NAME.net
    try:
        with open(join(NETS_DIR, net_name + '.net')) as f:
            lines = [line.strip() for line in f.read().rstrip().split('\n')]
    # some pseudo-nets not intended to be used on their own don't have the .net
    # extension
    except FileNotFoundError:
        with open(join(NETS_DIR, net_name)) as f:
            lines = [line.strip() for line in f.read().rstrip().split('\n')]
    isos = []

    # go line by line to find all explicitly mentioned isotopes
    in_isos = False
    for line in lines:
        # process lines within a call to add_isos or add_isos_and_reactions
        if in_isos:
            if line == ')':
                in_isos = False
            else:
                entries = line.rstrip(',').split()
                # ignore comments at the ends of lines
                if '!' in entries:
                    entries = entries[:entries.index('!')]
                # if single isotope, just add it.
                if len(entries) == 1:
                    isos.append(entries[0])
                # if a range of isotopes, add each one
                else:
                    element = entries[0]
                    a_start = int(entries[1])
                    a_end = int(entries[2])
                    for mass_number in range(a_start, a_end + 1):
                        isos.append(element + str(mass_number))
        # if we find a call to add_isos or add_isos_and_reactions,
        # start processing lines
        elif "add_isos" in line:
            if '(' in line and ')' in line:
                # extract comma-separated isotopes from inside the parentheses
                to_add = [iso.strip().replace("'", "").replace('"', '') for iso in line.split('(')[1].split(')')[0].split(',')]
                isos.extend(to_add)
                
            else:
                in_isos = True
        elif "add_iso" in line:
            if '(' in line and ')' in line:
                isos.append(line.split('(')[1].split(')')[0].split(',')[0])
        # if we find an include statement, recursively add isotopes from the
        # included net file
        elif line.startswith('include'):
            try:
                new_net = line.split()[1].replace("'", "").replace('"', '').replace('.net', '')
                isos.extend(isos_from_net(new_net))
            except IndexError as e:
                print("Could not parse include line in net file:")
                print(line)
                raise e

    # sort first by atomic number and then by mass number;
    # create a helper function to do this. Can't [easily] use a lambda because
    # the neutron is a special case
    def symbol_to_z_a(iso):
        """Take an isotope string and return a tuple of (Z, A)
        
        Z is atomic number (0 for neutrons) and A is mass number"""
        if iso == 'neut':
            return 0, 1
        return (symbol_to_z(re.split(r'(\d+)', iso)[0]), int(re.split(r'(\d+)', iso)[1]))
    # remove duplicates
    isos = list(set(isos))
    # sort in ascending (Z, A)
    isos.sort(key=symbol_to_z_a)
    return isos

if __name__ == "__main__":
    # if called as an executable, print the isotopes in the net file
    # provided as an argument
    net = argv[1]
    isos = isos_from_net(net)
    for i, iso in enumerate(isos):
        print(i + 1, iso)
