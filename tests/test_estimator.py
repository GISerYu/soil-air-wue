"""Mathematical unit tests; the small arrays below are synthetic, not observations."""
from pathlib import Path
import importlib.util
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('estimator', ROOT/'code/reproduce_models.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class EstimatorTests(unittest.TestCase):
    def test_bh_family(self):
        np.testing.assert_allclose(m.bh([.01,.04,.03,.20]), [.04,.053333333333,.053333333333,.20])

    def test_bh_order(self):
        p=np.array([.09,.01,.23,.001,.6])
        order=np.array([4,2,0,3,1])
        np.testing.assert_allclose(m.bh(p)[order], m.bh(p[order]))

    def test_shared_blocks(self):
        # Identical synthetic fits cancel exactly even with reversed block order.
        beta=np.zeros((11,3)); beta[7]=[.1,.03,.02]
        scores=np.zeros((3,11,3)); scores[:,7,:]=np.arange(9).reshape(3,3)/100
        a=dict(beta=beta,influence=scores,ids=np.array([1,2,3]),df=2)
        b=dict(beta=beta,influence=scores[::-1],ids=np.array([3,2,1]),df=2)
        out=m.estimate({'a':a,'b':b}, [('a',1),('b',-1)], 'WUE')
        self.assertAlmostEqual(out['beta'],0)
        self.assertAlmostEqual(out['se'],0)

    def test_cross_response_covariance(self):
        # Perfectly equal synthetic GPP and ET scores yield zero WUE uncertainty.
        beta=np.zeros((11,3)); scores=np.zeros((3,11,3))
        scores[:,7,0]=[1,2,3]; scores[:,7,2]=[1,2,3]
        f=dict(beta=beta,influence=scores,ids=np.array([1,2,3]),df=2)
        self.assertEqual(m.estimate({'f':f}, [('f',1)], 'WUE')['se'],0)

    def test_cluster_partition_point_estimate(self):
        # Real primary matrices: merge pairs of blocks without changing observations/design.
        path=ROOT/'data/model_inputs/primary_early.npz'
        original=m.fit(path)
        with np.load(path,allow_pickle=False) as z:
            n=len(z['ids'])
            groups=np.arange(n)//2
            unique=np.unique(groups)
            s=np.stack([z['s'][groups==g].sum(0) for g in unique])
            counts=np.stack([z['counts'][groups==g].sum(0) for g in unique])
        with tempfile.TemporaryDirectory() as folder:
            target=Path(folder)/'coarser.npz'
            np.savez_compressed(target,s=s,counts=counts,ids=unique)
            coarser=m.fit(target)
        np.testing.assert_allclose(original['beta'],coarser['beta'],atol=2e-11,rtol=2e-9)
        self.assertLess(coarser['n_clusters'],original['n_clusters'])
        self.assertEqual(coarser['n_cell_years'],original['n_cell_years'])

    def test_percent_conventions_differ(self):
        beta=.0123749012814278
        self.assertAlmostEqual(100*np.expm1(beta),1.245178719683156)
        self.assertNotAlmostEqual(100*beta,100*np.expm1(beta),places=3)


if __name__=='__main__':
    unittest.main()
