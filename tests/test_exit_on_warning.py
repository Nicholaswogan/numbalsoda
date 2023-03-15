import numpy as np
import numba as nb
import unittest

import numbalsoda


class TestNumbaLSODA(unittest.TestCase):

    def setUp(self):
        # numbalsoda
        @nb.cfunc(numbalsoda.lsoda_sig, boundscheck=False)
        def f_nb(t, u_, du_, p_):
            du_[0] = p_[0]*u_[0] * np.sin(t)

        self.f_nb = f_nb
        self.funcptr = self.f_nb.address

        self.atol = 1e-9
        self.rtol = 1e-9

    def test_exit_on_warning(self):
        """Test that setting exit_on_warning=True causes the solver to return
         prematurely terminate for an ill-conditioned ODE.

        """

        t_eval = np.linspace(100, 101, 11)
        u0 = np.array([1])
        args = np.array([1e15])

        usol_nb, success = numbalsoda.lsoda(self.funcptr, u0, t_eval, args,
                                            exit_on_warning=True,
                                            atol=self.atol, rtol=self.rtol)
        self.assertFalse(success)

    def test_solve_ode(self):
        """Test that the problem can be solved successfully if warnings
        are ignored.

        """

        t_eval = np.linspace(100, 101, 11)
        u0 = np.array([1])
        args = np.array([1e15])

        usol_nb, success = numbalsoda.lsoda(self.funcptr, u0, t_eval, args,
                                            atol=self.atol, rtol=self.rtol)
        self.assertTrue(success)
        return


if __name__ == "__main__":
    unittest.main()
