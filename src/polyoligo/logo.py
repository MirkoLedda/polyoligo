import sys

from . import _version


def main():
    argv = sys.argv[1:]

    print_logo = True
    if len(argv) > 0:
        if argv[0] == "cite":
            print("Ledda M., Cobo N., Lorant A., Hardigan M.A. and Knapp S.J., PolyOligo: A Bioinformatic Platform for "
                  "Identifying Target DNA Sequences for the Development of Sub-Genome Specific DNA Markers in "
                  "Polyploid/Complex Genomes. Poster presented at: Annual Conference of the American Society of "
                  "Horticultural Sciences; 2019 July 21-25; Las Vegas, NV, USA.\n")
            print_logo = False

    if print_logo:
        print(
            " ___     _       ___  _ _          \n"
            "| _ \___| |_  _ / _ \| (_)__ _ ___ \n"
            "|  _/ _ \ | || | (_) | | / _` / _ \ \n"
            "|_| \___/_|\_, |\___/|_|_\__, \___/ \n"
            "           |__/          |___/     v{}\n"
            "\n"
            "Design oligonucleotides for complex genomes. See https://github.com/MirkoLedda/polyoligo for more "
            "information.\n"
            "If you found this tool useful, please support us by citing -> see 'polyoligo cite'.\n"
            "Copyright 2019 to Mirko Ledda, under a BSD-2 License\n".format(_version.__version__)
        )


if __name__ == "__main__":
    main()
