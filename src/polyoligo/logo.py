import sys

from . import _version


def main():
    argv = sys.argv[1:]

    if len(argv) == 0:
        print(
            " ___     _       ___  _ _          \n"
            "| _ \___| |_  _ / _ \| (_)__ _ ___ \n"
            "|  _/ _ \ | || | (_) | | / _` / _ \ \n"
            "|_| \___/_|\_, |\___/|_|_\__, \___/ \n"
            "           |__/          |___/     v{}\n"
            "\n"
            "Design oligonucleotides for complex genomes. See https://github.com/MirkoLedda/polyoligo for more information.\n"
            "Please cite: TBA.\n"
            "(Bibtex reference available via 'polyoligo cite').\n"
            "\n"
            # "Available tools\n"
            # "   KASP       : polyoligo-kasp \n"
            # "   CAPS       : polyoligo-caps\n"
            # "   PCR/Sanger : polyoligo-pcr\n"
            # "   CRISPR/Cas9: polyoligo-crispr\n"
            # "\n"
            "Copyright 2019 to Mirko Ledda, under a BSD-2 License\n".format(_version.__version__)
        )
    else:
        print("Citation TBA.")


if __name__ == "__main__":
    main()
