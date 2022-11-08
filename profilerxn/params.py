import argparse


class Params:
    def __init__(self, **kwargs):
        coordsys_list = ["tric", "cart", "prim", "dlc", "hdlc"]
        self.fname = kwargs.get("xyzfile", None)
        self.charge = kwargs.get("charge", 0)
        self.mult = kwargs.get("mult", 1)
        self.method = kwargs.get("method", "b3lyp")
        self.basis = kwargs.get("basis", "6-31g*")
        self.frames = kwargs.get("frames", 20)
        self.coordsys = kwargs.get("coordsys", 'tric')
        self.springk = kwargs.get("springk", 0.1)
        self.workers = kwargs.get("workers", 2)
        self.equal_space = kwargs.get("equal_space", False)
        self.engine = kwargs.get("engine", "psi4")
        if self.coordsys not in coordsys_list:
            raise RuntimeError(
                "%s coordinate system is not available. Available coordinate systems: tric, cart, prim, dlc, hdlc"
                % self.coordsys
            )

        self.sp_spec = {
            "program": self.engine,
            "driver": "gradient",
            "method": self.method,
            "basis": self.basis,
        }

def str2bool(v):
    """ Allows command line options such as "yes" and "True" to be converted into Booleans. """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'on', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'off', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args(*args):
    parser = argparse.ArgumentParser(description="Parse arguments")
    parser.add_argument(
        "xyzfile",
        type=str,
        help="REQUIRED input xyzfile containing two optimized geometries",
    )
    parser.add_argument("--engine", type=str, help="QC software pacakge to perform energy/gradient calculation (default = psi4)")
    parser.add_argument("--charge", type=int, help="Molecular charge (default = 0)")
    parser.add_argument("--mult", type=int, help="Molecular multiplicity (default = 1)")
    parser.add_argument(
        "--method",
        type=str,
        help="Electronic structure method for singlepoint energy calculations (default = b3lyp)",
    )
    parser.add_argument(
        "--basis", type=str, help="Basis sets for the method (default = 6-31g*)"
    )
    parser.add_argument("--frames", type=int, help="Number of frames (default = 20)")
    parser.add_argument("--springk", type=float, help="Spring constant for NEB (default = 0.1)")
    parser.add_argument(
        "--coordsys",
        type=str,
        help="Coordinate system for interpolation (default = tric):\n"
        '"tric" for Translation-Rotation Internal Coordinates\n'
        '"cart" = Cartesian coordinate system\n'
        '"prim" = Primitive (a.k.a redundant internal coordinates)\n'
        '"dlc" = Delocalized Internal Coordinates,\n'
        '"hdlc" = Hybrid Delocalized Internal Coordinates\n',
    )
    parser.add_argument('--equal_space', type=str2bool,
                          help='Provide "yes" to space between frames equally at the end of the interpolation.\n')
    parser.add_argument(
        "--workers", type=int, help="Number of cores to be assigned to Psi4"
    )

    args_dict = {}
    for k, v in vars(parser.parse_args(*args)).items():
        if v is not None:
            args_dict[k] = v
    return args_dict
