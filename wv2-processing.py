#!/usr/bin/env python
# NOTE: use python3 in the shebang line above if targeting py3 and not py2-compatible
""" example main file with cmd line interface """

from argparse import ArgumentParser

from wv2_processing.main import main

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='short desc of projname goes here')

    parser.add_argument("source", help="directory to copy from")
    parser.add_argument('-l', '--log',
        help="desired filepath of log file",
        default="/var/opt/projectname/backup.log"
    )
    parser.add_argument('--rclonelog',
        help="desired path of rclone log file",
        default=None
    )
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="count",
                        default=0
    )

    args = parser.parse_args()

    main(args)
