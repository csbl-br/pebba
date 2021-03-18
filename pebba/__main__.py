from pebba.main import pebba


if __name__ == "__main__":
    import sys

    deg_file = sys.argv[1]
    gmt_file = sys.argv[2]

    pebba(
        deg_file,
        gmt_file,
    )
