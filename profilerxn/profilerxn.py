import os, sys, copy
from .params import Params, parse_args
import numpy as np
from geometric.molecule import Molecule
from geometric.molecule import EqualSpacing
from geometric.internal import (
    CartesianCoordinates,
    PrimitiveInternalCoordinates,
    DelocalizedInternalCoordinates,
)
from geometric.nifty import ang2bohr
from qcelemental.models import Molecule as qcel
from qcportal.neb import NEBKeywords
from qcfractal.snowflake import FractalSnowflake
import json


def check_record(record):
    status = record.status.split(".")[-1]
    if status != "complete":
        print(record)
        print(record.error)
        raise RuntimeError(
            "Failed task detected. See the printed error message above."
        )

def write_output(filename, dir, str):
    if not os.path.exists(dir):
        os.makedirs(dir)
    with open(os.path.join(dir, "%s" %filename), "w") as f:
        f.write(str)

class ProfileRXN:
    def __init__(self, params):

        file_dir = os.path.join(os.getcwd(), params.fname)
        if not os.path.exists(file_dir):
            raise RuntimeError("%s file can't be found." % file_dir)

        self.dir = os.getcwd()
        M = Molecule(file_dir)

        if M.elem == 0:
            raise RuntimeError(
                "xyz file format is not correct. Make sure there is no space between the two structures in the xyz file."
            )

        try:
            if len(M) > 2:
                print(
                    "The provided xyz file has more than two structures. The first and last geometries will be used."
                )
        except:
            raise RuntimeError(
                "Check your xyz file to make sure the two frames have the same numer of atoms"
            )

        self.elem = M.elem
        self.params = params

        self.s = FractalSnowflake(compute_workers=self.params.workers)
        self.client = self.s.client()

        self.sp_spec = params.sp_spec

        print(
            "Optimizing the input geometries with %s/%s"
            % (self.sp_spec["method"], self.sp_spec["basis"])
        )
        sp_spec = copy.copy(self.sp_spec)
        opt_spec = {"program": "geometric", "keywords": {}, "qc_specification": sp_spec}

        reac = qcel(
            symbols=self.elem,
            geometry=M.xyzs[0] * ang2bohr,
            molecular_charge=self.params.charge,
            molecular_multiplicity=self.params.mult,
        )
        prod = qcel(
            symbols=self.elem,
            geometry=M.xyzs[-1] * ang2bohr,
            molecular_charge=self.params.charge,
            molecular_multiplicity=self.params.mult,
        )

        _, ids = self.client.add_optimizations([reac, prod], **opt_spec)
        self.s.await_results()
        opt_recs = self.client.get_optimizations(ids)
        self.opt_dir = os.path.join(self.dir, "optimization_records")


        for i, rec in enumerate(opt_recs):
            check_record(rec)
            write_output("optimization_%i.out" %i, self.opt_dir, rec.stdout)
            #opt_output = open(os.path.join(self.opt_dir, "optimization_%i.out" %i), "w")
            #opt_output.write(rec.stdout)
            #opt_output.close()


        M.xyzs[0] = np.round(opt_recs[0].final_molecule.geometry / ang2bohr, 8)
        M.xyzs[-1] = np.round(opt_recs[-1].final_molecule.geometry / ang2bohr, 8)

        self.R_Energy = opt_recs[0].energies[-1]
        self.P_Energy = opt_recs[-1].energies[-1]
        if np.isclose(self.P_Energy, self.R_Energy, atol = 1e-06):
            print("WARNING: Optimized product and reactant have really close energies. It is possible that they converged to the same local minima.")
        # Coordinates in Angstrom
        M.build_topology()
        self.M = M
        self.interpolated_dict = {}
        self.Energy_dict = {}
        self.coordsys_dict = {
            "cart": (CartesianCoordinates, False, False),
            "prim": (PrimitiveInternalCoordinates, True, False),
            "dlc": (DelocalizedInternalCoordinates, True, False),
            "hdlc": (DelocalizedInternalCoordinates, False, True),
            "tric-p": (PrimitiveInternalCoordinates, False, False),
            "tric": (DelocalizedInternalCoordinates, False, False),
        }
        M_reac = self.M[0]
        M_prod = self.M[-1]
        M_reac.build_topology()
        M_prod.build_topology()

        if len(M_reac.molecules) >= len(M_reac.molecules):
            self.M_ini = copy.deepcopy(M_reac)
            self.M_fin = copy.deepcopy(M_prod)
        else:
            self.M_ini = copy.deepcopy(M_prod)
            self.M_fin = copy.deepcopy(M_reac)

        self.mix_xyz = self.M[0]
        self.fwd_M = self.M[0]
        self.bwd_M = self.M[0]

        self.reac = self.M_ini.xyzs[0].flatten() * ang2bohr
        self.prod = self.M_fin.xyzs[0].flatten() * ang2bohr

    def fill(self, cart1, cart2, ic, final_diff, mean_diff):
        reac_coords = cart1.copy()
        prod_coords = cart2.copy()
        filled_fwd_list = []
        filled_bwd_list = []
        M_ini = copy.deepcopy(self.M_ini)
        M_fin = copy.deepcopy(self.M_fin)
        M_ini.xyzs = [cart1.reshape(-1, 3) / ang2bohr]
        M_fin.xyzs = [cart2.reshape(-1, 3) / ang2bohr]
        CoordClass, connect, addacart = self.coordsys_dict[ic.lower()]
        nDiv = int(final_diff // mean_diff)
        nDiv += nDiv%2
        for i in range(nDiv // 2):
            IC_fwd = CoordClass(
                M_ini, build=True, connect=connect, addcart=addacart, constraints=None
            )
            IC_bwd = CoordClass(
                M_fin, build=True, connect=connect, addcart=addacart, constraints=None
            )
            if i == 0:
                dq_fwd = IC_fwd.calcDiff(reac_coords, prod_coords)
                dq_bwd = IC_bwd.calcDiff(prod_coords, reac_coords)
            else:
                dq_fwd = IC_fwd.calcDiff(new_bwd_coords, new_fwd_coords)
                dq_bwd = IC_bwd.calcDiff(new_fwd_coords, new_bwd_coords)
            interval = nDiv - 2 * i
            step_fwd = dq_fwd / interval
            step_bwd = dq_bwd / interval
            if interval <= 3:
                new_fwd_coords = IC_fwd.newCartesian(reac_coords, step_fwd)
                filled_fwd_list.append(new_fwd_coords)
                new_bwd_coords = IC_bwd.newCartesian(prod_coords, step_bwd)
                filled_bwd_list.append(new_bwd_coords)
                final_diff = np.linalg.norm(new_bwd_coords - new_fwd_coords)
                if final_diff / mean_diff > 1.0:
                    filled_list = self.fill(
                        new_fwd_coords, new_bwd_coords, ic, final_diff, mean_diff
                    )
                    filled_fwd_list += filled_list

                break
            else:
                new_fwd_coords = IC_fwd.newCartesian(reac_coords, step_fwd)
                new_bwd_coords = IC_bwd.newCartesian(prod_coords, step_bwd)
                reac_coords = new_fwd_coords.copy()
                prod_coords = new_bwd_coords.copy()

                filled_fwd_list.append(new_fwd_coords)
                filled_bwd_list.append(new_bwd_coords)

                M_ini.xyzs = [new_fwd_coords.reshape(-1, 3) / ang2bohr]
                M_fin.xyzs = [new_bwd_coords.reshape(-1, 3) / ang2bohr]
        filled_list = filled_fwd_list + filled_bwd_list[::-1]
        return filled_list

    def interpolate(self):
        ic = self.params.coordsys
        CoordClass, connect, addcart = self.coordsys_dict[ic.lower()]

        nDiv = self.params.frames
        nDiv += nDiv%2
        reac_coords = self.reac.copy()
        prod_coords = self.prod.copy()
        fwd_coord_list = [reac_coords]
        fwd_cart_diff = []
        bwd_coord_list = [prod_coords]
        bwd_cart_diff = []
        M_ini = copy.deepcopy(self.M_ini)
        M_fin = copy.deepcopy(self.M_fin)
        for i in range(nDiv // 2):
            IC_fwd = CoordClass(
                M_ini,
                build=True,
                connect=connect,
                addcart=addcart,
                constraints=None,
            )
            IC_bwd = CoordClass(
                M_fin,
                build=True,
                connect=connect,
                addcart=addcart,
                constraints=None,
            )

            if i == 0:
                dq_fwd = IC_fwd.calcDiff(self.prod, self.reac)
                dq_bwd = IC_bwd.calcDiff(self.reac, self.prod)
            else:
                dq_fwd = IC_fwd.calcDiff(new_bwd_coords, new_fwd_coords)
                dq_bwd = IC_bwd.calcDiff(new_fwd_coords, new_bwd_coords)

            interval = nDiv - 2 * i
            step_fwd = dq_fwd / interval
            step_bwd = dq_bwd / interval
            if interval == 2 :
                mean_diff = np.mean((fwd_diff_mean + bwd_diff_mean) / 2)
                new_fwd_coords = IC_fwd.newCartesian(reac_coords, step_fwd)
                fwd_coord_list.append(new_fwd_coords)
                new_bwd_coords = IC_bwd.newCartesian(prod_coords, step_bwd)
                bwd_coord_list.append(new_bwd_coords)
                final_diff = np.linalg.norm(fwd_coord_list[-1] - bwd_coord_list[-1])
                filled_list = self.fill(
                    new_fwd_coords, new_bwd_coords, ic, final_diff, mean_diff
                )
                fwd_coord_list += filled_list
                break
            else:
                new_fwd_coords = IC_fwd.newCartesian(reac_coords, step_fwd)
                new_bwd_coords = IC_bwd.newCartesian(prod_coords, step_bwd)

                fwd_diff = np.linalg.norm(
                    reac_coords.reshape(-1, 3) - new_fwd_coords.reshape(-1, 3),
                    axis=1,
                )
                fwd_cart_diff.append(fwd_diff)
                fwd_diff_mean = np.mean(fwd_cart_diff, axis=0)
                fwd_coord_list.append(new_fwd_coords)
                reac_coords = new_fwd_coords.copy()

                bwd_diff = np.linalg.norm(
                    prod_coords.reshape(-1, 3) - new_bwd_coords.reshape(-1, 3),
                    axis=1,
                )
                bwd_cart_diff.append(bwd_diff)
                bwd_diff_mean = np.mean(bwd_cart_diff, axis=0)
                bwd_coord_list.append(new_bwd_coords)
                prod_coords = new_bwd_coords.copy()

                M_ini.xyzs = [new_fwd_coords.reshape(-1, 3) / ang2bohr]
                M_fin.xyzs = [new_bwd_coords.reshape(-1, 3) / ang2bohr]

        coord_list = fwd_coord_list + bwd_coord_list[::-1]

        self.interpolated_dict["mix_" + ic] = np.array(coord_list)
        self.mix_xyz.xyzs = [coords.reshape(-1, 3) / ang2bohr for coords in coord_list]
        tot_len = len(self.mix_xyz)
        self.fwd_M.xyzs = [
            coords.reshape(-1, 3) / ang2bohr for coords in fwd_coord_list
        ]
        self.bwd_M.xyzs = [
            coords.reshape(-1, 3) / ang2bohr for coords in bwd_coord_list
        ]
        if self.params.equal_space:
            print("Spacing frames evenly.")
            final_M = EqualSpacing(self.mix_xyz)[
                np.array(
                    [
                        int(round(i))
                        for i in np.linspace(0, tot_len - 1, self.params.frames)
                    ]
                )
            ]
        else:
            final_M = self.mix_xyz[
                np.array(
                    [
                        int(round(i))
                        for i in np.linspace(0, tot_len - 1, self.params.frames)
                    ]
                )
            ]
        xyz_dir = os.path.join(self.dir, "interpolated")
        if not os.path.exists(xyz_dir):
            os.makedirs(xyz_dir)
        final_M.write(os.path.join(xyz_dir, "interpolated_%s.xyz" % ic))
        self.fwd_M.write(os.path.join(xyz_dir, "fwd_interpolated_%s.xyz" % ic))
        self.bwd_M.write(os.path.join(xyz_dir, "bwd_interpolated_%s.xyz" % ic))
        self.init_chain = final_M

    def run_NEB(self):
        chain = [
            qcel(
                symbols=M.elem,
                geometry=M.xyzs * ang2bohr,
                molecular_charge=self.params.charge,
                molecular_multiplicity=self.params.mult,
            )
            for M in self.init_chain
        ]
        neb_kw = NEBKeywords(
            images=self.params.frames,
            spring_constant=self.params.springk,
            optimize_ts=True,
            maximum_force=0.05,
            average_force=0.025,
            spring_type=0,
        )
        _, neb_id = self.client.add_nebs([chain], "geometric", self.sp_spec, neb_kw)
        self.s.await_results()

        neb_rec = self.client.get_nebs(neb_id[0])

        check_record(neb_rec)

        TS_output = open(os.path.join(self.opt_dir, "optimization_TS.out"), "w")
        TS_output.write(neb_rec.ts_optimization.stdout)
        TS_output.close()

        neb_dir = os.path.join(self.dir, "neb_result")
        if not os.path.exists(neb_dir):
            os.makedirs(neb_dir)
        neb_output = open(os.path.join(neb_dir, "neb.out"), "w")
        neb_output.write(neb_rec.stdout)
        neb_output.close()

        qcel_TS_M = neb_rec.ts_optimization.final_molecule
        self.TS_Energy = neb_rec.ts_optimization.energies[-1]

        TS_M = Molecule()
        TS_M.elem = qcel_TS_M.symbols.tolist()
        TS_M.comms = [
            "energy: %.14f calculated by %s" % (self.TS_Energy, self.params.engine)
        ]
        TS_M.xyzs = [np.round(qcel_TS_M.geometry / ang2bohr, 8)]

        self.TS_M = TS_M

        result_dir = os.path.join(self.dir, "final_result")
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)

        final_result_json = json.dumps(
            {
                "reactant": self.R_Energy,
                "transition": self.TS_Energy,
                "product": self.P_Energy,
            },
            indent=4,
        )

        write_output("final_energies.json", result_dir, final_result_json)
        #with open(os.path.join(result_dir, "final_energies.json"), "w") as f:
        #    f.write(final_result_json)

        self.M.comms = [
            "energy: %.14f calculated by %s" % (self.R_Energy, self.params.engine),
            "energy: %.14f calculated by %s" % (self.P_Energy, self.params.engine),
        ]
        self.M[0].write(os.path.join(result_dir, "reactant.xyz"))
        self.M[-1].write(os.path.join(result_dir, "product.xyz"))
        self.TS_M.write(os.path.join(result_dir, "transition.xyz"))


def main():
    args_dict = parse_args(sys.argv[1:])
    params = Params(**args_dict)

    PR = ProfileRXN(params)
    PR.interpolate()
    print("Interpolation is done.")
    print("NEB optimization is running now.")
    PR.run_NEB()
    print("Done!")


if __name__ == "__main__":
    main()
