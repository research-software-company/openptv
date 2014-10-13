"""
Check that the bindings do what they are expected to.
The test expects to be run from the py_bind/test/ directory for now,
using the nose test harness [1].

References:
[1] https://nose.readthedocs.org/en/latest/
"""

import unittest
import optv.imgcoord as ic
import numpy as np

class TestImgCoord(unittest.TestCase):

    @staticmethod
    def get_dummy_mm_np():
        ret = ic.pmm_np()
        ret.nlay = 3
        ret.n1 = 1.
        ret.n2 = np.array([1.49,0.,0])
        ret.d = np.array([5.,0.,0])
        ret.n3 = 1.33
        ret.lut = 1

        return ret

    @staticmethod
    def get_dummy_Interior():
        return {'xh':0. , 'yh':0. , 'cc':100.}


    @staticmethod
    def get_dummy_Glass():
        return {'vec_x': 0.0001,'vec_y': 0.00001,'vec_z': 1.}

    @staticmethod
    def get_dummy_Exterior():
        ret = ic.pExterior()
        ret.z0 = 100
        ret.dm = np.array([[1.0, 0.2, -0.3], 
            [0.2, 1.0, 0.0],
            [-0.3, 0.0, 1.0]])

        return ret 
        
    @staticmethod
    def get_dummy_addpar():
        return {'k1': 0.0, 'k2': 0.0, 'k3':0.0, 'p1':0.0, 'p2': 0.0, 'scx': 1.0, 'she': 0.0}

            
    @staticmethod
    def get_dummy_mmlut():
        mmlut = ic.pmmlut()
        mmlut.nr = 2
        mmlut.nz = 2
        mmlut.nw = 2
        mmlut.origin.x, mmlut.origin.y, mmlut.origin.z = 0.,0.,0.
        for i in range(mmlut.nr):
            for j in range(mmlut.nz):
                mmlut.data[i*mmlut.nz + j] = 1.
                
        return mmlut

    
    def test_img_coord(self):
        """Testing img_coord against data from test_imgcoord.c check testing."""
        
        mm = self.get_dummy_mm_np()        
        mmlut = self.get_dummy_mmlut()
        
        cal = ic.pCalibration()
        cal.ext_par = self.get_dummy_Exterior()
        cal.int_par = self.get_dummy_Interior()
        cal.glass_par = self.get_dummy_Glass()
        cal.added_par = self.get_dummy_addpar()


        
        imgcoord = ic.ImgCoord()
        imgcoord.force_set_calibration(cal)
        imgcoord.force_set_only_mm_np(mm)
        imgcoord.force_set_mmlut(mmlut)

        input_X = (100.0,100.0,0.0)
        int_cam = 0
        output_X = imgcoord.img_coord (input_X, int_cam)

        self.failUnlessAlmostEqual(
            np.max(
                np.abs(
                    np.array(output_X)-np.array((110.406944, 88.325788, 0.988076))
                )
            ),0. , places = 5)



