import pytest
import numpy as np 
import os 
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/..../')
from utils.io_couple import IO


class TestIO():

    def test_init(self,tmpdir):
        p = tmpdir.mkdir("test")
        io_obj = IO(
                    "kappa",
                    str(p),
                    str(p.join("test.txt")),
                    str(p.join("test2.txt")),
                    str(p.join("test3.txt")),
                    0,
                    True,
                    4)
        assert isinstance(io_obj._base_dir, str) 
        assert isinstance(io_obj._run_path, str) 
        assert isinstance(io_obj._run_name, str) 
        assert isinstance(io_obj.max_cycle, int) 
        assert isinstance(io_obj.cycle_counter, int) 

        assert io_obj._base_dir == p
        assert io_obj._run_path == os.path.join(p, "kappa")
        assert io_obj._k_src_dir == p.join("test.txt")
        assert io_obj._f_src_dir == p.join("test2.txt")
        assert io_obj._f_init_path == p.join("test3.txt")

    def test_create_dictionary(self, tmpdir):

        p = tmpdir.mkdir("test")
        io_obj = IO(
                    "kappa",
                    str(p),
                    str(p.join("test.txt")),
                    str(p.join("test2.txt")),
                    str(p.join("test3.txt")),
                    0,
                    4,
                    True)
        io_obj.createDirectoryOfOperation()
        io_obj.nextCyclePathManager()
        print(os.listdir(io_obj._run_path))
        print(io_obj.max_cycle)
        assert len(os.listdir(io_obj._run_path)) == 4
        assert len(os.listdir(io_obj.cycle_dump_path)) == 5


    def test_next_cycle_manager(self, tmpdir):

        p = tmpdir.mkdir("test")
        io_obj = IO(
                    "kappa",
                    str(p),
                    str(p.join("test.txt")),
                    str(p.join("test2.txt")),
                    str(p.join("test3.txt")),
                    0,
                    4,
                    True)

        io_obj.createDirectoryOfOperation()
        io_obj.nextCyclePathManager()
        assert len(os.listdir(io_obj.cycle_dump_path)) == 5
        assert io_obj.cycle_dump_path == os.path.join(p, 'kappa/CYCLE_0')
        assert io_obj.fluid_input_path == os.path.join(p, 'kappa/CYCLE_0/FLUID_INPUT/')
        assert io_obj.fluid_output_path == os.path.join(p, 'kappa/CYCLE_0/FLUID_OUTPUT/')
        assert io_obj.kinetic_input_path == os.path.join(p, 'kappa/CYCLE_0/KINETIC_INPUT/')
        assert io_obj.kinetic_output_path == os.path.join(p, 'kappa/CYCLE_0/KINETIC_OUTPUT/')
        assert io_obj.next_fluid_input_path == os.path.join(p, 'kappa/CYCLE_1/FLUID_INPUT/')

        io_obj.cycle_counter = 1
        io_obj.nextCyclePathManager()
        assert len(os.listdir(io_obj.cycle_dump_path)) == 5
        assert io_obj.cycle_dump_path == os.path.join(p, 'kappa/CYCLE_1')
        assert io_obj.fluid_input_path == os.path.join(p, 'kappa/CYCLE_1/FLUID_INPUT/')
        assert io_obj.fluid_output_path == os.path.join(p, 'kappa/CYCLE_1/FLUID_OUTPUT/')
        assert io_obj.kinetic_input_path == os.path.join(p, 'kappa/CYCLE_1/KINETIC_INPUT/')
        assert io_obj.kinetic_output_path == os.path.join(p, 'kappa/CYCLE_1/KINETIC_OUTPUT/')
        assert io_obj.next_fluid_input_path == os.path.join(p, 'kappa/CYCLE_2/FLUID_INPUT/')
        io_obj.cycle_counter = 2
        io_obj.nextCyclePathManager()
        assert len(os.listdir(io_obj.cycle_dump_path)) == 5
        assert io_obj.cycle_dump_path == os.path.join(p, 'kappa/CYCLE_2')
        assert io_obj.fluid_input_path == os.path.join(p, 'kappa/CYCLE_2/FLUID_INPUT/')
        assert io_obj.fluid_output_path == os.path.join(p, 'kappa/CYCLE_2/FLUID_OUTPUT/')
        assert io_obj.kinetic_input_path == os.path.join(p, 'kappa/CYCLE_2/KINETIC_INPUT/')
        assert io_obj.kinetic_output_path == os.path.join(p, 'kappa/CYCLE_2/KINETIC_OUTPUT/')
        assert io_obj.next_fluid_input_path == os.path.join(p, 'kappa/CYCLE_3/FLUID_INPUT/')
        io_obj.cycle_counter = 3
        io_obj.nextCyclePathManager()
        assert len(os.listdir(io_obj.cycle_dump_path)) == 5
        assert io_obj.cycle_dump_path == os.path.join(p, 'kappa/CYCLE_3')
        assert io_obj.fluid_input_path == os.path.join(p, 'kappa/CYCLE_3/FLUID_INPUT/')
        assert io_obj.fluid_output_path == os.path.join(p, 'kappa/CYCLE_3/FLUID_OUTPUT/')
        assert io_obj.kinetic_input_path == os.path.join(p, 'kappa/CYCLE_3/KINETIC_INPUT/')
        assert io_obj.kinetic_output_path == os.path.join(p, 'kappa/CYCLE_3/KINETIC_OUTPUT/')
        assert io_obj.next_fluid_input_path == os.path.join(p, 'kappa/CYCLE_4/FLUID_INPUT/')