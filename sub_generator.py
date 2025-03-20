####################################################################################
#                             Script by I. C. D. Merritt
#  Project: DE-TECT testing photochromes for novel photo-detectors
#  Function: Umbrella script
#  Created: 04/02/2024
#  Last edited: 04/02/2024
#  Usage: python
####################################################################################

import argparse
import os

from key_fns import list_substitution_library, run_input_checks, add_substitutions

script_location = os.path.realpath(os.path.dirname(__file__))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Substituent Generation Script')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-sub', '-S', action='store_true', help='List available files in substitution library')
    group.add_argument('-generate', '-G', type=str, help='Generate submission files for a set of substituted geometries')
    group.add_argument('-check_input', '-C', type=str, help='Check if given input file is acceptable.')
    args = parser.parse_args()

    if args.sub:
        list_substitution_library(script_location)
        exit()

    if args.generate:
        add_substitutions(args.generate, script_location)
        exit()

    if args.check_input:
        run_input_checks(args.check_input, script_location)
        exit()
