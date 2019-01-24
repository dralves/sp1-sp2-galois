import argparse
import pandas as pd;


parser = argparse.ArgumentParser(
    description='''Processes profiling outputs.''')

parser.add_argument(
    '-pf', '--profile_file',
    type=str,
    dest='profile_file',
    help='''The profile run output'''
)

parser.add_argument(
    '-s', '--sort_by',
    type=str,
    default='fixed_usec',
    dest='sort_by',
    help='''The column to sort the data by''')

parser.add_argument(
    '--sample_num',
    type=int,
    default=0,
    dest='sample_num',
    help='''Output every 'sample_num' value''')

if __name__ == "__main__":
    args = parser.parse_args()

    df = pd.read_csv(args.profile_file, sep=";")
    df = df.sort_values(by=args.sort_by)
    if args.sample_num != 0:
        df = df.iloc[::args.sample_num, :]

    df.to_csv(args.profile_file[0:-4] + "_sortby_" + args.sort_by + "_sample_n_" + str(args.sample_num) + ".csv")
