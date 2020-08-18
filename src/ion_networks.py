#!python

# local
try:
    from . import interface
except (ImportError, ModuleNotFoundError):
    import interface


def main():
    interface.CLI()


if __name__ == "__main__":
    main()
