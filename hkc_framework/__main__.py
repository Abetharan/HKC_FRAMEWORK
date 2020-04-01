from hkc_framework.coupler import Coupler
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Coupling of Hydrodynamic and Kinetic code. The two codes are specified in the input file.')
    parser.add_argument('-p', '--path', required=True, help = 'Give path to Input yml file.')
    args = vars(parser.parse_args())
    couple = Coupler(args['path'])
    couple.main()