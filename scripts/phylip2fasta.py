#!/usr/bin/env python3

import argparse
import sys
import os

parser = argparse.ArgumentParser(description="Converts newick trees in a single file to nexus format by simply just adding headings and terminal formatting.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree_in", dest="tree_in", required=True,
                                        help="Tree(s) file.")

requiredArgs.add_argument("-o", "--tree_out", dest="tree_out", required=True,
                                        help="Tree(s) file output.")

parser.add_argument("-r", "--overwrite", dest="remove", required=False, default=None, action="store_true",
                    help="Overwrite already existing output.")

args = parser.parse_args()

if os.path.exists(args.tree_out):
	print("\n  Warning! File", args.tree_out,"already exists.")
	if args.remove is not None:
		print("  Overwriting...\n")
		os.remove(args.tree_out)
	else:
		print("  Please choose other name for the output tree fle or consider using the option to overwrite (-r/--overwrite).\n")
		sys.exit(1)


with open(args.tree_out, "a") as outfile:
	print("#NEXUS", file=outfile)
	print("Begin trees;", file=outfile)
	c = 0
	for line in open(args.tree_in):
		c += 1
		print(f"\ttree tree_{c} = [&R] {line}", file=outfile, end="")
	print("end;", file=outfile)

