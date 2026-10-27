#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 2011.11.27 - Helgi I. Ingolfsson - Fix POPC and POPE (tail order is flipped in regular Martini 2.1 itp)
# 2015.07.01 - Alex H. de Vries - Adapt for lipidome (some lipid names changed, others not valid, 4 bead oleoyl tail)
# 2017.01.09 - Alex H. de Vries - Adapt for gromacs 5 and 2016

from math import sqrt
from os import remove, path
from sys import argv, stdout
import subprocess

def read_gro(file, atoms):
    """Parse a .gro file and return coordinates of two atoms."""
    line_counter = 0
    number_of_particles = 0
    first, second = [], []
    with open(file) as f:
        for line in f:
            if line_counter == 1:
                number_of_particles = int(line)
            elif 1 < line_counter < number_of_particles + 2:
                if line[10:15].strip() == atoms[0]:
                    first.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
                elif line[10:15].strip() == atoms[1]:
                    second.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
            line_counter += 1
    return [first, second]

if len(argv) != 11:
    print('''
Compute (second rank) order parameter, defined as:

  P2 = 0.5*(3*<cos²(theta)> - 1)

where "theta" is the angle between the bond and the bilayer normal.
P2 = 1      perfect alignment with the bilayer normal
P2 = -0.5   anti-alignment
P2 = 0      random orientation

Currently, the lipids DAPC, DLPC, DPPC, DSPC, DOPC, POPC, DLPE, DPPE, DSPE, DOPE, and POPE can be analyzed with this script.
It is quite simple to extend the range by editing this program. 
Search for DAPC and add similar lines for your lipid of interest. 

Usage: {0} <traj file> <tpr file> <initial time> <final time> <skip frames> <bilayer normal - xyz> <#lipids> <lipid type>

  Example:
    {0} my.xtc my.tpr 5000 10000 5 0 1 0 64 DPPC

Output: 
- P2 values in "order.dat"
- S-profile in "S-profile.dat"

Note: intermediate frames will be dumped as frame_dump_XXX.gro and deleted after processing.
'''.format(argv[0]))
    exit(0)

# === Input arguments ===
trajfile = argv[1]
tprfile = argv[2]
initial_time = int(argv[3])
final_time = int(argv[4])
traj_skip = int(argv[5])
orientation_of_bilayer_normal = [float(argv[6]), float(argv[7]), float(argv[8])]
number_of_lipids = int(argv[9])
lipid_type = argv[10]

# === Normalize bilayer normal vector ===
norm = sqrt(sum([x**2 for x in orientation_of_bilayer_normal]))
orientation_of_bilayer_normal = [x / norm for x in orientation_of_bilayer_normal]
stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" %
             tuple(orientation_of_bilayer_normal))

# === Define bond names per lipid ===
pc_bonds = " NC3-PO4 PO4-GL1 GL1-GL2 "
pe_bonds = " NH3-PO4 PO4-GL1 GL1-GL2 "
bond_dict = {
    "DAPC": pc_bonds + "GL1-D1A GL2-D1B D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B",
    "DLPC": pc_bonds + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B",
    "DOPC": pc_bonds + "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B",
    "DPPC": pc_bonds + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B",
    "DSPC": pc_bonds + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B",
    "POPC": pc_bonds + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B",
    "DLPE": pe_bonds + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B",
    "DOPE": pe_bonds + "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B",
    "DPPE": pe_bonds + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B",
    "DSPE": pe_bonds + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B",
    "POPE": pe_bonds + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B"
}
if lipid_type not in bond_dict:
    raise ValueError(f"Lipid type {lipid_type} not supported.")
bond_names = bond_dict[lipid_type]
output_legend = "  Frame" + bond_names

# === Write legend ===
stdout.write("\n " + output_legend + "\n")
stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n")

with open('order.dat', 'w') as output:
    output.write(output_legend + "\n")
    output.write("-"*(len(output_legend) - 1) + "\n")

    # === Extract frames ===
    stdout.write("Output all coordinate files\n")
    command = f"echo {lipid_type} | gmx trjconv -f {trajfile} -s {tprfile} -b {initial_time} -e {final_time} -sep -skip {traj_skip} -pbc whole -o frame_dump_.gro > /dev/null"
    print(command)
    subprocess.call(command, shell=True)

    # === Process each frame ===
    order_parameters = []
    file_count = 0
    bonds = [b.split("-") for b in bond_names.split()]

    stdout.write("Starting P2 calculation\n")
    while True:
        filename = f"frame_dump_{file_count}.gro"
        if not path.isfile(filename) or path.getsize(filename) == 0:
            break

        stdout.write(f"Processing snapshot {filename}\n")
        current_order_parameters = []

        for bond in bonds:
            first, second = read_gro(filename, bond)
            order_param = 0.0
            for i in range(number_of_lipids):
                vector = [first[i][j] - second[i][j] for j in range(3)]
                norm2 = sum([v**2 for v in vector])
                projection = sum([vector[j] * orientation_of_bilayer_normal[j] for j in range(3)])
                order_param += projection**2 / norm2
            p2 = 0.5 * (3.0 * (order_param / number_of_lipids) - 1.0)
            current_order_parameters.append(p2)

        order_parameters.append(current_order_parameters)

        # Write per-frame results
        results = f"{file_count:7d}" + "".join(f"{p:8.3f}" for p in current_order_parameters)
        stdout.write(" " + results + "\n")
        output.write(results + "\n")

        remove(filename)
        file_count += 1

    # === Final summary ===
    stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
    stdout.write("Snapshots analysis done.\n")
    stdout.write("Computing averages...\n")

    averaged_order_parameters = [sum(p[i] for p in order_parameters) / len(order_parameters)
                                 for i in range(len(bonds))]

    results = "average" + "".join(f"{p:8.3f}" for p in averaged_order_parameters)
    stdout.write("\n         " + bond_names + "\n")
    stdout.write("-"*(len(output_legend) - 1) + "\n")
    output.write("-"*(len(output_legend) - 1) + "\n")
    stdout.write(" " + results + "\n")
    output.write(results + "\n")
    stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")

    # === <Sn> for carbon chains only (from bond index 3 onward) ===
    ave_chain_s = sum(abs(val) for val in averaged_order_parameters[3:])
    sn_avg = ave_chain_s / (len(averaged_order_parameters) - 3)
    average_txt = f"Abs average order parameters for carbon chains <Sn> = {sn_avg:8.3f}\n\n"
    stdout.write(average_txt)
    output.write(average_txt)
    stdout.write('Results written in "order.dat".\n')

# === Write profile ===
with open('S-profile.dat', 'w') as profile:
    for i, val in enumerate(averaged_order_parameters):
        profile.write(f"{i+1:8.3f}   {val:8.3f}\n")
stdout.write('Order parameter profile written in "S-profile.dat".\n')

