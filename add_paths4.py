#!/usr/bin/python3

import os
import argparse as argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Transform a relative path to an absolute path')
    parser.add_argument('file', help='The file that needs its path defined.')
    parser.add_argument('directory', nargs='?', help='The current directory, relative to which the file may be found.', default=None)

    options = parser.parse_args()

    if options.directory is not None:
        os.chdir(options.directory)
    print(os.path.normpath(os.path.abspath(options.file)))
    
